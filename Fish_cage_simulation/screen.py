import numpy as np
import sys
from numpy import pi

np.set_printoptions(threshold=sys.maxsize)

row_air = 1.225  # [kg/m3]   공기 밀도
row_water = 1025.0  # [kg/m3]   물 밀도
gravity = 9.81

class ScreenModel:
    def __init__(self, mesh_elements, Sn, dw, rho,):
        self.two_elemnets = []
        self.dw = dw # 실제 그물발의 직경
        self.rho = rho # 재료밀도
        self.sn = Sn # 공극률 (잊지말고 meshinfo에 추가해야함)
        self.index = list(mesh_elements)

        # 0-based 원본 패널 저장
        converted_index = []
        for item in mesh_elements:
            # 중복 제거 + 1-based → 0-based
            uniq = [int(k) for k in dict.fromkeys(item)]  # 순서 유지 dedup
            converted_index.append([k-1 for k in uniq])
        self.original_elements_0b = converted_index  # e.g., [[0,1,2],[2,3,4,5]]

        for panel in self.original_elements_0b:
            if len(panel) <= 3:
                self.triangular_elements.append(panel)
                self._quad_of_sub.append(None)
            elif len(panel) == 4:
                a, b, c, d = panel
                # -1: 가상 중앙점 임시로 설정
                self.triangular_elements += [[a, b, -1], [b, c, -1], [c, d, -1], [d, a, -1]]
                self._quad_of_sub += [[a, b, c, d]] * 4

        self.hydro_dynamic_forces = np.zeros((len(self.triangular_elements), 3)) #힘 배열준비 [Fx1, Fy1, Fz1], [Fx2, Fy2, Fz2], ...
        self.hydro_z_forces = np.zeros((len(self.triangular_elements), 3))
        self.hydro_total_forces = np.zeors((len(self.triangular_elements), 3))

    def _resolve_triangle_points(self, tri, node_position, idx):
        """
        tri: [i,j,k] (k가 -1이면 중앙점)
        node_position: (N,3) 실노드 좌표
        idx: self.triangular_elements 내 현재 삼각형 인덱스
        return: p1, p2, p3 (모두 R^3)
        """
        i, j, k = tri
        if k != -1:  # 일반 삼각형
            return node_position[i], node_position[j], node_position[k]
        # 중앙점 삼각형: 부모 사각형의 코너 4개 평균
        a, b, c, d = self._quad_of_sub[idx]  # [a,b,c,d] 0-based
        center = (node_position[a] + node_position[b] +
                node_position[c] + node_position[d]) / 4.0
        return node_position[i], node_position[j], center

    def hydro_coefficients(self, point1, point2, point3, fluid_velocity, knot=False):
        a1 = point2 - point3
        a2 = point1 - point3
        unit_normal_vector = np.cross(a1, a2) / np.linalg.norm(np.cross(a1, a2))
        if np.dot(unit_normal_vector, fluid_velocity) < 0:
            unit_normal_vector = -unit_normal_vector
        surface_area = 0.5 * np.linalg.norm(np.cross(a1, a2))
        lift_vector = np.cross(np.cross(fluid_velocity, unit_normal_vector), fluid_velocity) / np.linalg.norm(
            np.cross(np.cross(fluid_velocity, unit_normal_vector), fluid_velocity) + np.finfo(np.float64).eps) # U는
        drag_vector = fluid_velocity / np.linalg.norm(fluid_velocity)
        
        coin_alpha = np.dot(unit_normal_vector, fluid_velocity) / np.linalg.norm(fluid_velocity) # 두 사이 각 구하기 ()
        inflow_angle = np.arccos(coin_alpha)        

        drag_coefficient, lift_coefficient = 0, 0 # 초기화 

        drag_coefficient = 0.04 + (
                -0.04 + 0.33 * self.sn + 6.54 * pow(self.sn, 2) - 4.88 * pow(self.sn, 3)) * np.cos(
            inflow_angle)
        lift_coefficient = (-0.05 * self.sn + 2.3 * pow(self.sn, 2) - 1.76 * pow(self.sn, 3)) * np.sin(
            2 * inflow_angle)
        
        return surface_area, drag_coefficient, lift_coefficient, drag_vector, lift_vector

    def force_on_element(self, node_position, velocity_fluid, velocity_structure=None, elevation=None):
        if velocity_structure is None:
            velocity_structure = np.zeros_like(node_position)

        num_element = len(self.triangular_elements)
        force_on_element = np.zeros((num_element, 3))
        eps = np.finfo(np.float64).eps

        # 🔹 (NEW) 사각형별 중앙점 속도 캐시: key = (a,b,c,d)  → value = v_center (R^3)
        quad_center_vel = {}

        for idx, tri in enumerate(self.triangular_elements):
            # 1) 좌표
            p1, p2, p3 = self._resolve_triangle_points(tri, node_position, idx)

            # 2) 상대속도
            u_f = velocity_fluid[idx]
            if tri[2] != -1:
                # 원래 삼각형: 세 노드 평균
                v_s = (velocity_structure[tri[0]] +
                    velocity_structure[tri[1]] +
                    velocity_structure[tri[2]]) / 3.0
            else:
                # 사각형 sub-tri: 중앙점 속도 = 4코너 평균 (캐시 활용)
                key = tuple(self._quad_of_sub[idx])   # (a,b,c,d)
                v_center = quad_center_vel.get(key)
                if v_center is None:
                    a, b, c, d = key
                    v_center = (velocity_structure[a] +
                                velocity_structure[b] +
                                velocity_structure[c] +
                                velocity_structure[d]) / 4.0
                    quad_center_vel[key] = v_center
                # 삼각형 꼭짓점(i, j, center)의 속도 평균
                v_s = (velocity_structure[tri[0]] +
                    velocity_structure[tri[1]] +
                    v_center) / 3.0

            rel = u_f - v_s

            area, Cd, Cl, drag_dir, lift_dir = self.hydro_coefficients(p1, p2, p3, rel)

            # 4) 밀도 선택
            rho = row_water
            if elevation is not None:
                elem_center = (p1 + p2 + p3) / 3.0
                if tri[2] != -1:
                    elev_vec = (elevation[tri[0]] + elevation[tri[1]] + elevation[tri[2]]) / 3.0
                else:
                    a, b, c, d = self._quad_of_sub[idx]
                    elev_center = (elevation[a] + elevation[b] + elevation[c] + elevation[d]) / 4.0
                    elev_vec = (elevation[tri[0]] + elevation[tri[1]] + elev_center) / 3.0
                elev_z = elev_vec[2] if np.ndim(elev_vec) >= 1 and len(np.atleast_1d(elev_vec)) >= 3 else float(elev_vec)
                rho = row_air if elem_center[2] > elev_z else row_water

            # 5) 힘 크기
            U = np.linalg.norm(rel) + eps
            FD_mag = 0.5 * rho * Cd * area * (U**2)
            FL_mag = 0.5 * rho * Cl * area * (U**2)

            # 6) 요소 힘
            FD = FD_mag * drag_dir
            FL = FL_mag * lift_dir
            force_on_element[idx] = (FD + FL) / 2.0  # Cheng 스케일 유지

        # 7) 상태 업데이트
        self.hydro_dynamic_forces = force_on_element
        self.hydro_total_forces   = force_on_element
        return force_on_element
        

    def cal_buoy_force(self, node_position, elevation):
        """
        요소별 정수압(부력) 힘을 계산해 self.hydro_z_forces에 저장하고 반환.
        - 요소 중심 z < 수면 z 일 때만 작용 (아니면 0)
        - 압력 p = rho_water * g * depth  (depth = elev_z - z_center)
        - 힘 = p * area * n_hat  (n_hat의 z성분이 음수면 뒤집어 위(+z) 성분이 되도록)
        """
        num_element = len(self.triangular_elements)
        buoy_forces = np.zeros((num_element, 3), dtype=float)

        if elevation is None:
            # 수면 정보 없으면 부력 0 (필요 시 상수 수면을 넘기세요)
            self.hydro_z_forces = buoy_forces
            # total은 동유체력만 유지
            self.hydro_total_forces = getattr(self, "hydro_dynamic_forces", np.zeros_like(buoy_forces))
            return buoy_forces

        for idx, tri in enumerate(self.triangular_elements):
            # 요소 좌표
            p1, p2, p3 = self._resolve_triangle_points(tri, node_position, idx)
            # 요소 면적과 법선(속도와 무관하게 기하로 계산)
            a1 = p2 - p3
            a2 = p1 - p3
            cross = np.cross(a1, a2)
            area = 0.5 * (np.linalg.norm(cross) + np.finfo(np.float64).eps)
            n_hat = cross / (2.0 * area)  # 단위법선 (cross/(||cross||))
            # 위쪽(+z)으로 작용하도록 법선 정방향 선택
            if n_hat[2] < 0.0:
                n_hat = -n_hat

            # 요소 중심 z
            elem_center = (p1 + p2 + p3) / 3.0

            # 요소 주변 수면고(z) 산정 (동유체력 블록과 동일 규칙)
            if tri[2] != -1:
                elev_vec = (elevation[tri[0]] + elevation[tri[1]] + elevation[tri[2]]) / 3.0
            else:
                a, b, c, d = self._quad_of_sub[idx]
                elev_center = (elevation[a] + elevation[b] + elevation[c] + elevation[d]) / 4.0
                elev_vec = (elevation[tri[0]] + elevation[tri[1]] + elev_center) / 3.0

            elev_z = elev_vec[2] if (np.ndim(elev_vec) >= 1 and len(np.atleast_1d(elev_vec)) >= 3) else float(elev_vec)

            # 수면 아래만 정수압 작용
            depth = elev_z - float(elem_center[2])
            if depth > 0.0 and area > 0.0:
                p = row_water * gravity * depth  # 정수압
                buoy_forces[idx] = p * area * n_hat  # 요소 부력 벡터

        self.hydro_z_forces = buoy_forces

        # 총합력 업데이트(동유체력 + 정수압)
        dyn = getattr(self, "hydro_dynamic_forces", np.zeros_like(buoy_forces))
        self.hydro_total_forces = dyn + buoy_forces
        return buoy_forces


    def distribute_force(self, number_of_node):
        """
        요소 힘(self.hydro_total_forces가 있으면 그걸, 없으면 동유체력)을
        노드 힘 (Nnode x 3)으로 분배해서 반환.
        - 일반 삼각형: 각 꼭짓점에 1/3씩
        - sub-tri([i,j,-1]): 실제 노드 i,j에 1/2씩 (가상 중앙점은 분배하지 않음)
        """
        N = int(number_of_node)
        node_forces = np.zeros((N, 3), dtype=float)

        # 분배 대상 선택: 총합력이 있으면 그걸 우선
        elem_forces = getattr(self, "hydro_total_forces", None)
        if elem_forces is None or elem_forces.shape[0] != len(self.triangular_elements):
            elem_forces = getattr(self, "hydro_dynamic_forces", np.zeros((len(self.triangular_elements), 3)))

        for idx, tri in enumerate(self.triangular_elements):
            f = elem_forces[idx]
            i, j, k = tri
            if k != -1:
                # 3개 꼭짓점 균등 분배
                share = f / 3.0
                if 0 <= i < N: node_forces[i] += share
                if 0 <= j < N: node_forces[j] += share
                if 0 <= k < N: node_forces[k] += share
            else:
                # sub-tri: 실제 두 노드에 1/2씩
                share = f / 2.0
                if 0 <= i < N: node_forces[i] += share
                if 0 <= j < N: node_forces[j] += share
                # 중앙점은 실제 노드가 아니므로 분배하지 않음

        return node_forces



if __name__ == "__main__":
    pass

