import numpy as np
import sys
from numpy import pi

np.set_printoptions(threshold=sys.maxsize)

row_air = 1.225  # [kg/m3]   공기 밀도
row_water = 1025.0  # [kg/m3]   물 밀도
gravity = 9.81

class ScreenModel:
    def __init__(self, mesh_elements, Sn, dw, rho, lines_netting = None):
        self.two_elemnets = []
        self.dw = dw # 실제 그물발의 직경
        self.rho = rho # 재료밀도
        self.sn = Sn # 공극률 (잊지말고 meshinfo에 추가해야함)
        self.index = list(mesh_elements)

        self.triangular_elements = []
        self._quad_of_sub = []

        # 0-based 원본 패널 저장
        converted_index = []
        for item in mesh_elements:
            # 중복 제거 + 1-based → 0-based
            uniq = [int(k) for k in dict.fromkeys(item)]  # 순서 유지 dedup
            converted_index.append([k-1 for k in uniq])
        self.original_elements_0b = converted_index  # e.g., [[0,1,2],[2,3,4,5]]

        # ✅ Lines_netting을 이용해 4코너 순서 자동 정렬
        self._edges = None
        if lines_netting:
            es = set()
            for e in lines_netting:
                if len(e) != 2:
                    continue
                i, j = e[0] - 1, e[1] - 1
                if i == j:
                    continue
                es.add((i, j) if i < j else (j, i))
            self._edges = es

        # ✅ 기존 패널 분해 로직 (단, 순서 정렬 추가)
        for panel in self.original_elements_0b:
            if len(panel) <= 3:
                self.triangular_elements.append(panel)
                self._quad_of_sub.append(None)
            elif len(panel) == 4:
                ordered = self._order_quad_with_edges(panel)
                if ordered is None:
                    ordered = panel  # fallback (Lines_netting이 없을 경우)
                a, b, c, d = ordered
                self.triangular_elements += [
                    [a, b, -1],
                    [b, c, -1],
                    [c, d, -1],
                    [d, a, -1],
                ]
                self._quad_of_sub += [[a, b, c, d]] * 4

        self.hydro_dynamic_forces = np.zeros((len(self.triangular_elements), 3)) #힘 배열준비 [Fx1, Fy1, Fz1], [Fx2, Fy2, Fz2], ...
        self.hydro_z_forces = np.zeros((len(self.triangular_elements), 3))
        self.hydro_total_forces = np.zeros((len(self.triangular_elements), 3))

    # -----------------------------------------------------------
    #  Lines_netting 기반 정렬 함수 
    # -----------------------------------------------------------
    def _order_quad_with_edges(self, nodes4):
        if self._edges is None:
            return None

        # 각 노드의 연결 정보
        nbr = {u: set() for u in nodes4}
        for i in range(4):
            for j in range(i + 1, 4):
                u, v = nodes4[i], nodes4[j]
                key = (u, v) if u < v else (v, u)
                if key in self._edges:
                    nbr[u].add(v)
                    nbr[v].add(u)

        # 연결 불완전 시 실패
        if any(len(nbr[u]) != 2 for u in nodes4):
            return None

        # 둘레 순회 (연결 따라 돌기)
        start = nodes4[0]
        order = [start]
        prev, cur = None, start
        for _ in range(3):
            nxt = [n for n in nbr[cur] if n != prev][0]
            order.append(nxt)
            prev, cur = cur, nxt
        return order

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

        # 사각형별 중앙점 속도 캐시: key = (a,b,c,d)  → value = v_center (R^3)
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
            force_on_element[idx] = (FD + FL)

        # 7) 상태 업데이트
        self.hydro_dynamic_forces = force_on_element
        self.hydro_total_forces   = force_on_element
        return force_on_element
        

    def cal_buoy_force(self, node_position, elevation):
        """
        elevation: np.array shape (N,) or (N,3) 모두 지원
        """
        num_element = len(self.triangular_elements)
        force_on_element = np.zeros((num_element, 3))
        eps = np.finfo(np.float64).eps

        def _elev_z(elev, idx):
            """elevation[idx]에서 z 성분/스칼라 높이 추출"""
            v = elev[int(idx)]
            v = np.asarray(v)
            if v.ndim == 0:
                return float(v)
            return float(v[2] if v.size >= 3 else v.item())

        def _submerge_ratio(zvals):
            # 물속(>0)=1, 수면(=0)=0.5, 공기(<0)=0
            zvals = np.asarray(zvals, dtype=float)
            r = ((zvals > 0).astype(float) + 0.5 * (np.abs(zvals) <= eps).astype(float)).mean()
            return max(0.0, min(1.0, float(r)))

        for index, tri in enumerate(self.triangular_elements):
            # 1) 좌표
            p1, p2, p3 = self._resolve_triangle_points(tri, node_position, index)

            # 2) 면적 및 체적
            element_area = self.hydro_coefficients(p1, p2, p3, np.array([1, 0, 0]))[0]
            element_volume = element_area * self.sn * self.dw * 0.25 * pi

            # 3) 수면고 z 추출 (중앙점이면 4 코너 평균)
            i, j, k = tri
            if k != -1:
                e0 = _elev_z(elevation, i)
                e1 = _elev_z(elevation, j)
                e2 = _elev_z(elevation, k)
            else:
                a, b, c, d = self._quad_of_sub[index]
                e_center = (_elev_z(elevation, a) + _elev_z(elevation, b) +
                            _elev_z(elevation, c) + _elev_z(elevation, d)) / 4.0
                e0 = _elev_z(elevation, i)
                e1 = _elev_z(elevation, j)
                e2 = e_center

            # 4) 잠김 비율 평가 (elev_z - node_z)
            z_rel = [e0 - p1[2], e1 - p2[2], e2 - p3[2]]
            ratio_water = _submerge_ratio(z_rel)

            # 5) 유체 밀도 보간 및 부력
            rho_fluid = row_air * (1 - ratio_water) + row_water * ratio_water
            Fz = element_volume * gravity * (rho_fluid)
            force_on_element[index] = np.array([0.0, 0.0, Fz])

        self.hydro_static_forces = force_on_element
        return force_on_element


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
                # 원래 삼각형인 그물코
                share = f / 3.0
                if 0 <= i < N: node_forces[i] += share
                if 0 <= j < N: node_forces[j] += share
                if 0 <= k < N: node_forces[k] += share
            else:
                # 4등분한 그물코의 두 노드에 힘을 분배
                share = f / 2.0
                if 0 <= i < N: node_forces[i] += share
                if 0 <= j < N: node_forces[j] += share
                # 중앙점은 실제 노드가 아니므로 분배하지 않음

        return node_forces

if __name__ == "__main__":

    node_position = np.array([
        [0.0, 1.0, 0.0],  # Node 0
        [0.0, 0.0, 0.0],  # Node 1
        [0.7071, 1.0, -0.7071],  # Node 2
        [0.7071, 0.0, -0.7071],  # Node 3
    ])

    lines_netting = [[1,2],[2,4],[4,3],[3,1]]  # 1-based, 테두리 엣지

    surfs_netting = [[1, 2, 3, 4]]  # 1-based

    screen = ScreenModel(surfs_netting, Sn=0.3, dw=0.002, rho=1140.0, lines_netting=lines_netting)

    velocity_fluid = np.tile(np.array([1.0, 0.0, 0.0]), (len(screen.triangular_elements), 1))

    # ✅ 구조속도는 0으로 가정
    velocity_structure = np.zeros_like(node_position)

    # ✅ 수면고 (z=0.5m)
    elevation = np.full((len(node_position), 3), [0, 0, 0.5])

    # ✅ 부력 계산
    f_buoy = screen.cal_buoy_force(node_position, elevation)
    print("\n[부력 per element]\n", f_buoy)

    # ✅ 동유체력 계산
    f_dyn = screen.force_on_element(node_position, velocity_fluid, velocity_structure, elevation)
    print("\n[동유체력(Drag+Lift) per element]\n", f_dyn)