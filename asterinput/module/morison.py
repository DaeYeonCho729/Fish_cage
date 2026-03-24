import numpy as np
import sys
from numpy import pi
import os
np.set_printoptions(threshold=sys.maxsize)

row_air = 1.225  # [kg/m3]   공기 밀도
row_water = 1025.0  # [kg/m3]   물 밀도
gravity = 9.81

class MorisonModel:
    def __init__(self, line_elements, dw, rho, th, node_ids=None):

        converted = []
        for e in (line_elements or []):
            if e is None or len(e) < 2:
                continue
            i = int(e[0]) - 1
            j = int(e[1]) - 1
            if i == j:
                continue
            converted.append((i, j))

        self.line_elements = converted

        self.python_to_med = None
        if node_ids is not None:
            self.python_to_med = node_ids[:]  # copy 저장
        self.onedimension_elements = []
        self._quad_of_sub = []  
        self.two_elemnets = []
        self.dw = dw # 외부 지름
        self.rho = rho # 재료밀도
        self.dr = dw - 2 * th # 내부직경
        Ne = len(self.line_elements)
        self.hydro_dynamic_forces = np.zeros((Ne, 3), dtype=float)
        self.hydro_static_forces = np.zeros((Ne, 3), dtype=float)
        self.hydro_total_forces = np.zeros((Ne, 3), dtype=float)

    def output_hydro_element(self):
        """
        확인용 코드
        """
        return self.line_elements

    def hydro_coefficients(self, point1, point2, fluid_velocity, tri_nodes=None):
        eps = np.finfo(np.float64).eps
        tiny = 1e-12

        a = point2 - point1
        L = np.linalg.norm(a)
        U = np.linalg.norm(fluid_velocity)

        unit_normal_vector = a / (L)
        if np.dot(unit_normal_vector, fluid_velocity) < 0:
            unit_normal_vector = -unit_normal_vector

        # drag/lift 방향
        drag_vector = fluid_velocity / (U)
        lift_vector = np.cross(np.cross(fluid_velocity, unit_normal_vector), fluid_velocity) / np.linalg.norm(
            np.cross(np.cross(fluid_velocity, unit_normal_vector), fluid_velocity) + eps) 
        lift_vector = -lift_vector

        # 각도
        alpha = np.dot(unit_normal_vector, fluid_velocity) / (U)
        alpha = float(np.clip(alpha, -1.0, 1.0))
        attack_angle_radian = np.arccos(alpha)
        attack_angle_deg = np.degrees(attack_angle_radian)
        drag_coefficient, lift_coefficient = 0, 0 # 초기화 

        # Cd
        if 0 <= attack_angle_deg <= 10:
            drag_coefficient = 0.00185 * attack_angle_deg + 0.033352
        elif 10 < attack_angle_deg <= 20:
            drag_coefficient = 0.005986 * attack_angle_deg - 0.01081
        elif 20 < attack_angle_deg <= 40:
            drag_coefficient = 0.015606 * attack_angle_deg - 0.2101
        elif 40 < attack_angle_deg <= 60:
            drag_coefficient = 0.0239 * attack_angle_deg - 0.54333
        elif 60 < attack_angle_deg <= 70:
            drag_coefficient = 0.011907 * attack_angle_deg + 0.183297
        elif 70 < attack_angle_deg <= 90:
            drag_coefficient = 0.008774 * attack_angle_deg + 0.402622

        # Cl
        if 0 <= attack_angle_deg <= 5:
            lift_coefficient = 0.001831 * attack_angle_deg + 0.020675
        elif 5 < attack_angle_deg <= 10:
            lift_coefficient = 0.004817 * attack_angle_deg + 0.004268
        elif 10 < attack_angle_deg <= 20:
            lift_coefficient = 0.008873 * attack_angle_deg - 0.03891
        elif 20 < attack_angle_deg <= 30:
            lift_coefficient = 0.008873 * attack_angle_deg - 0.11956
        elif 30 < attack_angle_deg <= 43:
            lift_coefficient = 0.013325 * attack_angle_deg - 0.13140
        elif 43 < attack_angle_deg <= 50:
            lift_coefficient = 0.008639 * attack_angle_deg + 0.070035
        elif 50 < attack_angle_deg <= 55:
            lift_coefficient = 0.002428 * attack_angle_deg + 0.383586
        elif 55 < attack_angle_deg <= 60:
            lift_coefficient = -0.00503 * attack_angle_deg + 0.797889
        elif 60 < attack_angle_deg <= 65:
            lift_coefficient = -0.01468 * attack_angle_deg + 1.38167
        elif 65 < attack_angle_deg <= 72:
            lift_coefficient = -0.02828 * attack_angle_deg + 2.274163
        elif 72 < attack_angle_deg <= 85:
            lift_coefficient = -0.01127 * attack_angle_deg + 1.046885
        elif attack_angle_deg > 85:
            lift_coefficient = -0.00703 * attack_angle_deg + 0.689754

        return float(L), float(drag_coefficient), float(lift_coefficient), drag_vector, lift_vector

    def force_on_element(self, node_position, velocity_fluid, velocity_structure=None, elevation=None):
        if velocity_structure is None:
            velocity_structure = np.zeros_like(node_position)

        num_element = len(self.line_elements)
        force_on_element = np.zeros((num_element, 3), dtype=float)
        eps = np.finfo(np.float64).eps

        vfluid = np.asarray(velocity_fluid, dtype=float)
        pos = np.asarray(node_position, dtype=float)
        vel_s = np.asarray(velocity_structure, dtype=float)

        N = int(pos.shape[0])
        tiny = 1e-12

        def _elev_z(elev, idx):
            if elev is None:
                return 0.0
            v = elev[int(idx)]
            v = np.asarray(v)
            if v.ndim == 0:
                return float(v)
            return float(v[2] if v.size >= 3 else v.item())

        def _segment_submerge_ratio(p1, p2, e1, e2):
            """
            선요소 잠김 비율(0~1)
            z_rel = elev - z
            z_rel > 0 : 수면 아래
            z_rel < 0 : 수면 위
            """
            z1 = float(e1 - p1[2])
            z2 = float(e2 - p2[2])

            # 완전 수중
            if z1 >= 0.0 and z2 >= 0.0:
                return 1.0
            # 완전 수상
            if z1 <= 0.0 and z2 <= 0.0:
                return 0.0

            # 부분 잠김: 선형 보간
            denom = (z1 - z2)
            if abs(denom) < tiny:
                return 0.5

            t = z1 / denom  # z_rel=0 되는 위치
            t = max(0.0, min(1.0, t))

            # p(t)=p1+t(p2-p1), t=0 -> p1, t=1 -> p2
            # z1>0, z2<0 이면 p1 쪽이 물속
            if z1 > 0.0 and z2 < 0.0:
                return t
            # z1<0, z2>0 이면 p2 쪽이 물속
            elif z1 < 0.0 and z2 > 0.0:
                return 1.0 - t
            else:
                return 0.0

        for idx, (n1, n2) in enumerate(self.line_elements):
            if n1 < 0 or n2 < 0 or n1 >= N or n2 >= N:
                force_on_element[idx] = 0.0
                continue

            p1 = pos[n1]
            p2 = pos[n2]

            # 유체 속도(요소 대표값)
            if vfluid.shape[0] == N:
                u_f = 0.5 * (vfluid[n1] + vfluid[n2])
            else:
                u_f = vfluid[idx]

            # 잠김 비율
            if elevation is None:
                ratio_water = 1.0
            else:
                e1 = _elev_z(elevation, n1)
                e2 = _elev_z(elevation, n2)
                ratio_water = _segment_submerge_ratio(p1, p2, e1, e2)

            # 수면 위면 유체속도 0, 부분잠김이면 비율만큼 감소
            u_f = ratio_water * u_f

            # 구조 속도(요소 대표값)
            v_s = 0.5 * (vel_s[n1] + vel_s[n2])

            while np.linalg.norm(v_s) > np.linalg.norm(u_f) and np.linalg.norm(u_f) > tiny:
                v_s *= 0.1

            if np.dot(v_s, u_f) > 0:
                rel = u_f - v_s
            else:
                rel = u_f

            U = float(np.linalg.norm(rel))
            if U < tiny or ratio_water <= 0.0:
                force_on_element[idx] = 0.0
                continue

            # 계수 + 방향벡터
            L, Cd, Cl, drag_dir, lift_dir = self.hydro_coefficients(p1, p2, rel, tri_nodes=(n1, n2))

            # 유체 밀도도 잠김 비율로 보간
            rho = row_air * (1.0 - ratio_water) + row_water * ratio_water

            # projected area
            A = float(self.dw) * float(L)

            FD = 0.5 * rho * Cd * A * (U ** 2) * drag_dir
            FL = 0.5 * rho * Cl * A * (U ** 2) * lift_dir

            force_on_element[idx] = FD + FL

        self.hydro_dynamic_forces = force_on_element
        self.hydro_total_forces = force_on_element
        return force_on_element

    def cal_buoy_force(self, node_position, elevation):
        """
        선요소(로프/플로터/바텀링) 부력(정수력) 계산
        - elevation: shape (N,) 또는 (N,3) 지원
        - 잠김 비율 기반으로 공기/물 밀도 보간
        """
        pos = np.asarray(node_position, dtype=float)
        Nnode = pos.shape[0]
        Ne = len(self.line_elements)
        force_on_element = np.zeros((Ne, 3), dtype=float)
        tiny = 1e-12

        dws = float(self.dw)

        def _elev_z(elev, idx):
            v = elev[int(idx)]
            v = np.asarray(v)
            if v.ndim == 0:
                return float(v)
            return float(v[2] if v.size >= 3 else v.item())

        def _segment_submerge_ratio(p1, p2, e1, e2):
            z1 = float(e1 - p1[2])
            z2 = float(e2 - p2[2])

            if z1 >= 0.0 and z2 >= 0.0:
                return 1.0
            if z1 <= 0.0 and z2 <= 0.0:
                return 0.0

            denom = (z1 - z2)
            if abs(denom) < tiny:
                return 0.5

            t = z1 / denom
            t = max(0.0, min(1.0, t))

            if z1 > 0.0 and z2 < 0.0:
                return t
            elif z1 < 0.0 and z2 > 0.0:
                return 1.0 - t
            else:
                return 0.0

        for eidx, (i, j) in enumerate(self.line_elements):
            if not (0 <= i < Nnode and 0 <= j < Nnode):
                continue

            p1 = pos[i]
            p2 = pos[j]
            a = p2 - p1
            L = float(np.linalg.norm(a))
            if L < tiny:
                continue

            V_g = 0.25 * pi * (dws**2 - self.dr**2) * L

            V_b = 0.25 * pi * (dws**2) * L

            e1 = _elev_z(elevation, i)
            e2 = _elev_z(elevation, j)
            ratio_water = _segment_submerge_ratio(p1, p2, e1, e2)

            rho_fluid = row_air * (1.0 - ratio_water) + row_water * ratio_water

            buoy = V_b * gravity * rho_fluid
            weight = V_g * gravity * self.rho

            Fz = buoy - weight
            force_on_element[eidx] = np.array([0.0, 0.0, Fz], dtype=float)

        self.hydro_static_forces = force_on_element
        return force_on_element

    def distribute_force(self, number_of_node):
        N = int(number_of_node)
        node_forces = np.zeros((N, 3), dtype=float)

        elem_forces = getattr(self, "hydro_total_forces", None)
        if elem_forces is None or elem_forces.shape[0] != len(self.line_elements):
            elem_forces = getattr(self, "hydro_dynamic_forces",
                                np.zeros((len(self.line_elements), 3), dtype=float))

        for idx, (i, j) in enumerate(self.line_elements):
            f = elem_forces[idx]
            share = 0.5 * f

            if 0 <= i < N:
                node_forces[i] += share
            if 0 <= j < N:
                node_forces[j] += share

        return node_forces


if __name__ == "__main__":
    pass