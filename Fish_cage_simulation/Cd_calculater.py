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
if __name__ == "__main__":
    pass