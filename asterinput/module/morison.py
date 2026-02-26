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

        vfluid = np.asarray(velocity_fluid, float)

        # -------------------- [B옵션 추가: 1-based 자동 감지/변환] --------------------
        N = int(np.asarray(node_position).shape[0])  # node_position의 노드 개수
        one_based = False
        if num_element > 0:
            # self.line_elements에서 최대 인덱스를 보고 판단
            try:
                max_idx = max((n1 if n1 > n2 else n2) for (n1, n2) in self.line_elements)
                if max_idx >= N:
                    one_based = True
            except Exception:
                one_based = False
        # ---------------------------------------------------------------------------

        for idx, (n1, n2) in enumerate(self.line_elements):
            # -------------------- [B옵션 추가: 변환 + 범위검사] --------------------
            if one_based:
                n1 -= 1
                n2 -= 1

            if n1 < 0 or n2 < 0 or n1 >= N or n2 >= N:
                force_on_element[idx] = 0.0
                continue
            # ----------------------------------------------------------------------
            p1 = node_position[n1]
            p2 = node_position[n2]
            tiny = 1e-12
            # 유체 속도(요소 대표값)
            if vfluid.shape[0] == node_position.shape[0]:
                u_f = 0.5 * (vfluid[n1] + vfluid[n2])
            else:
                u_f = vfluid[idx]

            # 구조 속도(요소 대표값)
            v_s = 0.5 * (velocity_structure[n1] + velocity_structure[n2])

            # Screen에서 했던 안전장치
            while np.linalg.norm(v_s) > np.linalg.norm(u_f):
                v_s *= 0.1
            if np.dot(v_s, u_f) > 0:
                rel = u_f - v_s
            else :
                rel = u_f - v_s * 0


            # 계수(각도기반 Cd/Cl) + 방향벡터
            L, Cd, Cl, drag_dir, lift_dir = self.hydro_coefficients(p1, p2, rel, tri_nodes=(n1, n2))

            # 밀도(유체)
            rho = row_water

            # 원통 projected area
            A = float(self.dw) * float(L)

            U = float(np.linalg.norm(rel))

            FD = 0.5 * rho * Cd * A * (U ** 2) * drag_dir
            FL = 0.5 * rho * Cl * A * (U ** 2) * lift_dir

            force_on_element[idx] = FD + FL

        self.hydro_dynamic_forces = force_on_element
        self.hydro_total_forces = force_on_element
        return force_on_element

    def cal_buoy_force(self, node_position, elevation):
        """
        선요소(로프/플로터/바텀링) 부력(정수력) 계산
        - element volume: V = (pi/4) * d^2 * L
        - elevation: shape (N,) 또는 (N,3) 지원
        - 잠김판정: 두 끝 노드의 (elev_z - node_z) 부호로 ratio_water 결정(0~1)
        """
        pos = np.asarray(node_position, dtype=float)
        Nnode = pos.shape[0]
        Ne = len(self.line_elements)
        force_on_element = np.zeros((Ne, 3), dtype=float)
        tiny = 1e-12

        dws = self.dw

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

            V_b = 0.25 * pi * (dws ** 2) * L

            rho_fluid = row_water

            buoy = V_b * gravity * (rho_fluid)
            weight = V_g * gravity * (self.rho)
            
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