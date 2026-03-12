import numpy as np
import sys
from numpy import pi
import os
np.set_printoptions(threshold=sys.maxsize)

row_air = 1.225  # [kg/m3]   공기 밀도
row_water = 1025.0  # [kg/m3]   물 밀도
gravity = 9.81

class ScreenModel:
    def __init__(self,rr, mesh_elements, Sn, dw, rho, lines_netting = None, wake_origin = None, node_ids=None):

        self.python_to_med = None
        if node_ids is not None:
            self.python_to_med = node_ids[:]  # copy 저장
        self.triangular_elements = []
        self._quad_of_sub = []  
        self.rr = rr
        self.two_elemnets = []
        self.dw = dw # 실제 그물발의 직경
        self.dr = dw * rr # 보간 직경
        self.rho = rho # 재료밀도
        self.sn = Sn # 공극률 (잊지말고 meshinfo에 추가해야함)
        # mesh_elements 처리하여 original_elements_0b와 태그 분리
        converted_index = []
        self._net_tags = []
        for item in mesh_elements:
            if item and isinstance(item[-1], str):
                tag = item[-1]             # 태그 분리
                nodes = item[:-1]
            else:
                tag = None
                nodes = item
            # 중복 제거 및 1-베이스→0-베이스 변환
            uniq_nodes = [int(k) for k in dict.fromkeys(nodes)]
            converted_index.append([k - 1 for k in uniq_nodes])
            self._net_tags.append(tag)    # 각 패널의 태그 저장 (태그 없으면 None)
        self.original_elements_0b = converted_index
        self._tag_map = {frozenset(panel): (self._net_tags[i] if i < len(self._net_tags) else None)
                         for i, panel in enumerate(self.original_elements_0b)}
        
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

        # --- wake effect를 위한 변수 ---
        self.wake_initialized = False
        self.wake_element_indexes = []  # wake 판정된 sub-triangle 인덱스들
        # (n,1)로 두면 브로드캐스팅이 안전함
        self.wake_reduction_factors = np.ones((len(self.triangular_elements), 1), dtype=float)
        self._flow_dir = np.array([1.0, 0.0, 0.0], dtype=float)
        self._flow_speed = 0.0
        self._wake_origin = wake_origin

    #  Lines_netting 기반 정렬 함수 
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

    def output_hydro_element(self):
        """
        확인용 코드
        """
        return self.triangular_elements

    def hydro_coefficients(self, point1, point2, point3, fluid_velocity,
                        tri_nodes=None):

        # 좌표
        eps = np.finfo(np.float64).eps
        tiny = 1e-12

        # cross 계산
        a1 = point2 - point3
        a2 = point1 - point3
        normal = np.cross(a1, a2)
        norm_n = np.linalg.norm(normal)

        # =======================================================
        # 🔥 문제 발생시(= cross 길이 0) → 로그만 기록하고 끝
        # =======================================================
        if norm_n < tiny:
            with open("hydro_error_log.txt", "a", encoding="utf-8") as f:
                f.write("\n==== CROSS ZERO ERROR ====\n")

                # python 인덱스 출력
                if tri_nodes is not None:
                    f.write(f"triangle python indices : {tri_nodes}\n")

                    # MED 번호 출력 가능하면 출력
                    if getattr(self, "python_to_med", None) is not None:
                        try:
                            med_ids = [self.python_to_med[i] for i in tri_nodes if i >= 0]
                            f.write(f"triangle MED nodes    : {med_ids}\n")
                        except:
                            f.write("triangle MED nodes    : mapping failed\n")

                f.write(f"a1 : {a1.tolist()}\n")
                f.write(f"a2 : {a2.tolist()}\n")
                f.write(f"cross(a1,a2) : {normal.tolist()}\n")
                f.write(f"norm(cross)  : {norm_n}\n")
                f.write("==========================\n")

            # 방어 x → 그냥 None 반환 (네가 원하는 대로)
            return 0.0, 0.0, 0.0, np.zeros(3), np.zeros(3)

        Urel = float(np.linalg.norm(fluid_velocity))
        unit_normal_vector = normal / (norm_n)
        if np.dot(unit_normal_vector, fluid_velocity) < 0:
            unit_normal_vector = -unit_normal_vector
        surface_area = 0.5 * np.linalg.norm(np.cross(a1, a2))
        
        lift_vector = np.cross(np.cross(fluid_velocity, unit_normal_vector), fluid_velocity) / np.linalg.norm(
            np.cross(np.cross(fluid_velocity, unit_normal_vector), fluid_velocity) + eps) # U는
        drag_vector = fluid_velocity / np.linalg.norm(fluid_velocity)
        
        attack_angle_deg = np.dot(unit_normal_vector, fluid_velocity) / np.linalg.norm(fluid_velocity) # 두 사이 각 구하기 ()
        attack_angle = np.arccos(attack_angle_deg)   #라디안    
        inflow_angle = (pi/2) - attack_angle #라디안

        drag_coefficient, lift_coefficient = 0, 0 # 초기화 

        # -----------------------------
        # 3) 해수 20도 동점성계수
        # -----------------------------
        nu = 1.05e-6   # [m^2/s]

        # Eq. (11)
        Re = (self.dw * Urel) / (nu * (1.0 - self.sn) + eps)

        x = np.log10(Re)

        # Eq. (10): C_D^(circ:cyl)
        cd_circ_cyl = (
            -78.46675
            + 254.73873 * x
            - 327.88640 * (x ** 2)
            + 223.64577 * (x ** 3)
            - 87.92234 * (x ** 4)
            + 20.00769 * (x ** 5)
            - 2.44894 * (x ** 6)
            + 0.12479 * (x ** 7)
        )

        # -----------------------------
        # 4) cd 계산: Section 2.5.1
        #    cd = CN(0), Eq. (9) with y=0
        # -----------------------------
        cd = (
            cd_circ_cyl * self.sn * (2.0 - self.sn)
            / (((1.0 - self.sn) ** 2) *2)
        )

        # -----------------------------
        # 5) cl 계산: Section 2.5.1
        #    CN(pi/4)=0.5*cd
        #    CT(pi/4) from Eq. (7)
        #    cl from Eq. (2)
        # -----------------------------
        cn_45 = 0.5 * cd
        y_45 = 0.25 * np.pi

        ct_45 = y_45 * ((4.0 * cn_45) / (8.0 + cn_45))

        cl = (cn_45 - ct_45) / (2**0.5)

        drag_coefficient = cd * np.cos(inflow_angle)
        lift_coefficient = cl * np.sin(inflow_angle * 2)

        return (
            float(surface_area),
            float(drag_coefficient),
            float(lift_coefficient),
            drag_vector,
            lift_vector,
        )

    def force_on_element(self, node_position, velocity_fluid, velocity_structure=None, elevation=None):
        if velocity_structure is None:
            velocity_structure = np.zeros_like(node_position)

        num_element = len(self.triangular_elements)
        force_on_element = np.zeros((num_element, 3))
        eps = np.finfo(np.float64).eps

        # 사각형별 중앙점 속도 저장
        quad_center_vel = {}

        for idx, tri in enumerate(self.triangular_elements):
            # 1) 좌표
            p1, p2, p3 = self._resolve_triangle_points(tri, node_position, idx)

            # 2) 상대속도
            u_f = velocity_fluid[idx]

            if self.wake_initialized:
                r = float(self.wake_reduction_factors[idx])   # (n,1) → 스칼라
                u_f = r * u_f

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

            while np.linalg.norm(v_s) > np.linalg.norm(u_f):
                v_s *= 0.1
            if np.dot(v_s, u_f) > 0:
                rel = u_f - v_s
            else :
                rel = u_f - v_s * 0

            area, Cd, Cl, drag_dir, lift_dir = self.hydro_coefficients(p1, p2, p3, rel, tri_nodes=tri)

            # 4) 밀도 선택
            rho = row_water

            # 5) 힘 크기
            U = np.linalg.norm(rel)
            FD_mag = 0.5 * rho * Cd * area * (U**2) * drag_dir
            FL_mag = 0.5 * rho * Cl * area * (U**2)

            # 6) 요소 힘
            FD = FD_mag
            FL = FL_mag * lift_dir
            force_on_element[idx] = (FD + FL)
            
            #print(f"[Elem {idx:02d}] Cd={Cd:.4f}, Cl={Cl:.4f}, |U|={U:.3f}, Area={area:.5f}")

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
            element_area = self.hydro_coefficients(p1, p2, p3, np.array([1, 0, 0]), tri_nodes=tri)[0]
            tag = None
            if self._net_tags:  # 태그 정보가 있을 경우에만 처리
                if tri[2] == -1:
                    # 사각 패널에서 생성된 sub-triangle인 경우, 원본 4코너 셋으로 태그 결정
                    tag = self._tag_map.get(frozenset(self._quad_of_sub[index]), None)
                else:
                    # 원본이 삼각 패널인 경우, 해당 삼각형 노드 셋으로 태그 결정
                    tag = self._tag_map.get(frozenset(tri), None)
                    
            # 태그에 따라 부력 계산에 사용할 직경 선택 ('b'이면 dw, 그 외는 dr)
            if tag == 'b':
                element_volume = element_area * self.sn * self.dw * 0.25 * pi
            else:
                element_volume = element_area * self.sn * self.dr * 0.25 * pi

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
            Fz = element_volume * gravity * (rho_fluid - self.rho)
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

#########################################################################
#                           Wake effect                                 #
#########################################################################
    def init_wake(self, initial_node_position, flow_vec):
        """
        초기 형상과 유동벡터로 wake 판정만 고정해두는 초기화.
        - initial_node_position : (N,3) 초기 노드 좌표
        - flow_vec              : (3,) 유동 벡터 (절대속도)
        - mode                  : 추후 확장 대비 문자열 태그
        """
        pos = np.asarray(initial_node_position, dtype=float)
        flow_vec = np.asarray(flow_vec, dtype=float)

        self._flow_speed = float(np.linalg.norm(flow_vec))
        self._flow_dir = flow_vec / (self._flow_speed + np.finfo(np.float64).eps)
        self._wake_origin = pos.mean(axis=0)

        # 초기 형상 기준: (패널 중심 - origin) · flow_dir > 0 → wake 영역
        idxs = []
        for i, tri in enumerate(self.triangular_elements):
            # 내부 보조함수로 안전하게 좌표 해석
            p1, p2, p3 = self._resolve_triangle_points(tri, pos, i)
            center = (p1 + p2 + p3) / 3.0
            if np.dot(center - self._wake_origin, self._flow_dir) > 0.0:
                idxs.append(i)

        self.wake_element_indexes = idxs
        self.wake_reduction_factors[:] = 1.0  # 아직 감속율 r은 모두 1.0
        self.wake_initialized = True

    def panel_cd(self, idx, node_position, flow_velocity_vec, velocity_structure=None):
        """
        idx: self.triangular_elements 기준 삼각 서브패널 인덱스
        flow_velocity_vec: (3,) 이 패널에 적용될 유체 속도 벡터
        반환: float(c_d)
        """
        if velocity_structure is None:
            velocity_structure = np.zeros_like(node_position)

        tri = self.triangular_elements[idx]
        # 좌표 해석(중앙 가상노드 케이스 포함) — 기존 보조함수 사용
        p1, p2, p3 = self._resolve_triangle_points(tri, node_position, idx)

        # 요소 구조 속도(네 현재 규칙 유지: 실노드 3개 평균 or 가상중심 포함 평균)
        if tri[2] != -1:
            v_s = (velocity_structure[tri[0]] +
                   velocity_structure[tri[1]] +
                   velocity_structure[tri[2]]) / 3.0
        else:
            # 사각 -> 삼각 분해 케이스: 중앙점은 가상노드 → 사각형 4개 실노드 평균을 중심속도로 사용
            a, b, c, d = self._quad_of_sub[idx]
            v_center = (velocity_structure[a] + velocity_structure[b] +
                        velocity_structure[c] + velocity_structure[d]) / 4.0
            v_s = (velocity_structure[tri[0]] +
                   velocity_structure[tri[1]] +
                   v_center) / 3.0

        rel_u = np.asarray(flow_velocity_vec, dtype=float) - v_s
        _, c_d, _, _, _ = self.hydro_coefficients(p1, p2, p3, rel_u, tri_nodes=tri)
        return float(c_d)

    def _r_from_cd(self, c_d):
        """
        c_d -> wake 감속율 r(0~1)로 변환. 
        """
        r = 1.0 - 0.46 * float(c_d)
        return max(0.0, min(1.0, r))
    
    def update_wake_reduction(self, node_position, velocity_structure=None):
        """
        init_wake()로 지정된 wake 요소들만 r 갱신.
        아직 force_on_element에는 적용하지 않는다(다음 단계).
        """
        if not self.wake_initialized:
            return
        if velocity_structure is None:
            velocity_structure = np.zeros_like(node_position)

        flow_vec = self._flow_dir * self._flow_speed  # 초기 유동을 유지

        for i in self.wake_element_indexes:
            c_d = self.panel_cd(i, node_position, flow_vec, velocity_structure)
            self.wake_reduction_factors[i] = self._r_from_cd(c_d)

if __name__ == "__main__":
    pass