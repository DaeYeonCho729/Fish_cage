import sys
import json
import numpy as np
import os
import math, json, glob

# ✅ 누락된 위젯/코어 import 보강
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QFrame, QWidget,
    QVBoxLayout, QHBoxLayout, QPushButton, QLineEdit, QLabel,
    QFormLayout, QTabWidget, QDialog, QComboBox, QListWidget,
    QComboBox, QListWidget, QRadioButton, QButtonGroup, QMessageBox, QCheckBox, QSplitter, QMenu, QAction)


from PyQt5.QtCore import Qt, QTimer, QEvent, QSize
from PyQt5.QtGui import QDoubleValidator

# ✅ VTK import를 환경별로 안전하게
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor

try:
    from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
except Exception:
    import vtk
    vtkInteractorStyleTrackballCamera = vtk.vtkInteractorStyleTrackballCamera

try:
    from vtkmodules.vtkInteractionWidgets import vtkOrientationMarkerWidget
    from vtkmodules.vtkRenderingAnnotation import vtkAxesActor
except Exception:
    import vtk
    vtkOrientationMarkerWidget = vtk.vtkOrientationMarkerWidget
    vtkAxesActor = vtk.vtkAxesActor

import vtk

# vtkUnsignedCharArray도 환경별 안전하게
try:
    from vtkmodules.vtkCommonCore import vtkUnsignedCharArray as _VtkUCA
except Exception:
    _VtkUCA = vtk.vtkUnsignedCharArray
vtkUnsignedCharArray = _VtkUCA

from geometry_generator import generate_cage_frame
from net_generator import generate_net_nodes_and_edges
from mooring_generator import generate_mooring_frame
from Floating_Collar import generate_outer_floater
from export_vtu import export_all_to_vtu

BOOT_FLAG = -1
LINKS_DIR = os.path.join("Fish_Cage", "links")
def ensure_links_dir():
    os.makedirs(LINKS_DIR, exist_ok=True)

#카메라 설정
class PanOnlyStyle(vtkInteractorStyleTrackballCamera):
    def Rotate(self):
        pass 

class MooringRigWindow(QDialog):
    """
    Mooring Frame의 선택 노드 기준으로
    - Buoy tether(상향)
    - Distance rope(수평, 방사 방향)
    - Buoy line(수면 방향 상향)
    - Anchor line(하향)
    을 한 번에 생성하여 links JSON으로 저장
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Mooring Rigging")
        self.resize(960, 640)

        # ---- UI (split: left VTK, right form) ----
        root = QHBoxLayout(self)

        # Left: VTK preview
        self.vtkWidget = QVTKRenderWindowInteractor(self)
        self.renderer = vtk.vtkRenderer()
        self.vtkWidget.GetRenderWindow().AddRenderer(self.renderer)
        self.renderer.SetBackground(0.18, 0.18, 0.18)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.iren.Initialize()
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        root.addWidget(self.vtkWidget, stretch=2)
        
        # Right: controls
        panel = QWidget(self)
        form = QFormLayout(panel)

        # 대상 mooring 노드
        self.node_combo = QComboBox()
        form.addRow(QLabel("Mooring node"), self.node_combo)

        # 파라미터
        self.in_buoy_up = QLineEdit("3.0")
        self.in_dist_len = QLineEdit("10.0")
        self.in_buoyline_up = QLineEdit("2.0")
        self.in_anchor_down = QLineEdit("15.0")
        self.in_azimuth_deg = QLineEdit("")
        self.anchor_incl_deg = QLineEdit("45")
        self.anchor_incl_deg.setValidator(QDoubleValidator(0.0, 89.9, 1))
        self.anchor_incl_deg.setToolTip("수직 기준 각도 (0°=수직, 90°=수평)")

        form.addRow(QLabel("Buoy 1 (m)"), self.in_buoy_up)
        row_b2 = QWidget(); hb_b2 = QHBoxLayout(row_b2); hb_b2.setContentsMargins(0,0,0,0)
        self.cb_buoy2 = QCheckBox("use")
        self.cb_buoy2.setChecked(False)
        hb_b2.addWidget(self.in_buoyline_up); hb_b2.addWidget(self.cb_buoy2)
        form.addRow(QLabel("Buoy 2 (m)"), row_b2)

        form.addRow(QLabel("Distance rope length (m)"), self.in_dist_len)
        form.addRow(QLabel("Anchor-line down (m)"), self.in_anchor_down)
        form.addRow(QLabel("Anchor Inclination (deg)"), self.anchor_incl_deg)

        self.in_anchor_a = QLineEdit("1") 
        self.in_anchor_b = QLineEdit("1")
        w_ab = QWidget(); lay_ab = QHBoxLayout(w_ab); lay_ab.setContentsMargins(0,0,0,0)
        lay_ab.addWidget(self.in_anchor_a); lay_ab.addWidget(QLabel(":")); lay_ab.addWidget(self.in_anchor_b)
        form.addRow(QLabel("Anchor split (A:B)"), w_ab)

        form.addRow(QLabel("Azimuth (deg, optional)"), self.in_azimuth_deg)

        # ▼ UniversalLinkWindow처럼 Add / List / Delete
        self.btn_add = QPushButton("Add")
        self.sel_list = QListWidget()
        self.btn_del = QPushButton("Delete")
        form.addRow(self.btn_add)
        form.addRow(self.sel_list)
        form.addRow(self.btn_del)

        # Save/Close
        btns = QHBoxLayout()
        self.btn_save = QPushButton("Save")
        self.btn_close = QPushButton("Close")
        btns.addWidget(self.btn_save); btns.addWidget(self.btn_close)
        form.addRow(btns)

        root.addWidget(panel, stretch=1)

        # 데이터 로드
        self.sel_nodes = []           # 전역 gidx 목록 (Add로 누적)
        self.display_points = []      # 미리보기 좌표(로컬)
        self.display_map = {}         # 글로벌 gidx -> 로컬 인덱스
        self.items = []  # ← 각 Add 시점의 {gidx, buoy_up, dist_len, buoyline_up, anchor_down, azimuth_deg} 누적

        self._last_gidx = None  # ✅ 선택이 비어도 미리보기 유지용 저장소
        self.load_data()
        self.populate_nodes()
        self.build_preview_data()     # ⬅️ 미리보기 포인트 채우기
        self.display_vtk()            # ⬅️ 초기 렌더
        self.reload_existing_rig()   # ★ 추가

        # 시그널
        self.node_combo.currentIndexChanged.connect(self.display_vtk)
        self.btn_add.clicked.connect(self.add_current_node)
        self.btn_del.clicked.connect(self.del_selected_node)
        self.btn_save.clicked.connect(self.save_rig)
        self.btn_close.clicked.connect(self.close)
        # Buoy2 enable/disable
        def _apply_buoy2_enabled(on: bool):
            self.in_dist_len.setEnabled(on)
            self.in_buoyline_up.setEnabled(on)
        self._apply_buoy2_enabled = _apply_buoy2_enabled
        if hasattr(self, "cb_buoy2"):
            self.cb_buoy2.stateChanged.connect(lambda *_: (self._apply_buoy2_enabled(self.cb_buoy2.isChecked()), self.display_vtk()))
            self._apply_buoy2_enabled(self.cb_buoy2.isChecked())

        for w in (self.in_buoy_up, self.in_dist_len, self.in_buoyline_up, self.in_anchor_down, self.anchor_incl_deg, self.in_azimuth_deg):
            w.editingFinished.connect(self.display_vtk)
            
    def _get_anchor_incl_deg(self):
        try:
            return max(0.0, min(89.9, float(self.anchor_incl_deg.text())))
        except:
            return 0.0

    # ✅ 현재 콤보 선택이 없으면 마지막 선택을 반환
    def _current_gidx(self):
        g_raw = self.node_combo.currentData()
        try:
            g = int(str(g_raw)) if g_raw is not None else None
        except Exception:
            g = None
        if g is not None:
            self._last_gidx = g
        return self._last_gidx


    def reload_existing_rig(self):
        """
        links/MOORING_RIG.json을 읽어 self.items를 복구한다.
        - 선호: ANCHOR_LINE_A + ANCHOR_LINE_B (A/B 분할)
        - 폴백: ANCHOR_LINE (구버전 단일)
        - BUOY_TETHER, DISTANCE_ROPE, BUOY_LINE은 필수
        """
        # ----- 경로 -----
        try:
            links_dir = LINKS_DIR  # 프로젝트 전역 상수 있으면 사용
        except NameError:
            links_dir = os.path.join(os.path.dirname(__file__), "links")
        path = os.path.join(links_dir, "MOORING_RIG.json")
        if not os.path.isfile(path):
            return  # 파일 없으면 조용히 종료

        # ----- 로드 -----
        with open(path, "r", encoding="utf-8") as f:
            data = json.load(f)

        baseN      = int(data.get("baseN", 0))
        extra_pts  = data.get("extra_points", [])
        polylines  = data.get("polylines", [])
        if not polylines:
            return

        # ----- 유틸 -----
        def coord_of(node_id: int):
            """전역 node_id -> (x,y,z)"""
            if node_id < baseN:
                return tuple(self._get_point_by_gidx(int(node_id)))
            j = node_id - baseN
            if 0 <= j < len(extra_pts):
                p = extra_pts[j]
                return (float(p[0]), float(p[1]), float(p[2]))
            raise KeyError(f"node {node_id} coord not found")

        def seglen3D(p, q):
            return ((p[0]-q[0])**2 + (p[1]-q[1])**2 + (p[2]-q[2])**2) ** 0.5

        def best_ratio(la, lb, max_den=12):
            """
            A/B 실수 비(la:lb)를 1..max_den 정수로 근사.
            la+lb가 0이면 (1,1). 한쪽이 극단적으로 짧아도 1은 보장.
            """
            total = la + lb
            if total <= 1e-12:
                return 1, 1
            fa = la / total
            best = (1, max_den)
            best_err = abs(fa - best[0] / (best[0] + best[1]))
            for a in range(1, max_den+1):
                for b in range(1, max_den+1):
                    err = abs(fa - a / (a + b))
                    if err < best_err - 1e-9:
                        best, best_err = (a, b), err
            return best

        # 그룹 색인
        by_group = {}
        for pl in polylines:
            by_group.setdefault(pl.get("group", ""), []).append(pl)

        # BUOY_TETHER: (base gidx, buoy idx) 매핑
        #   표준: [gidx, i_buoy], 역순 저장 대비
        tether_of_buoy = {}  # buoy_idx -> base_gidx
        for pl in by_group.get("BUOY_TETHER", []):
            n0, n1 = map(int, pl["nodes"])
            if n0 < baseN and n1 >= baseN:
                tether_of_buoy[n1] = n0
            elif n1 < baseN and n0 >= baseN:
                tether_of_buoy[n0] = n1

        # 1) Buoy2 있는 형식: DISTANCE_ROPE 기준으로 한 세트씩 복구
        restored = []
        for dr in by_group.get("DISTANCE_ROPE", []):
            d0, d1 = map(int, dr["nodes"])
            # ✅ 새 포맷: [base gidx(<baseN), i_dist_end(>=baseN)]
            gidx, i_dist, i_buoy = None, None, None
            if (d0 < baseN) ^ (d1 < baseN):
                gidx  = d0 if d0 < baseN else d1
                i_dist = d1 if d0 < baseN else d0
                # base에 매단 부이 찾기(테더 역참조)
                for pl in by_group.get("BUOY_TETHER", []):
                    n0, n1 = map(int, pl["nodes"])
                    if n0 == gidx and n1 >= baseN:
                        i_buoy = n1; break
                    if n1 == gidx and n0 >= baseN:
                        i_buoy = n0; break
                if i_buoy is None:
                    continue
            else:
                # ⬅️ 구포맷: [i_buoy(>=baseN), i_dist_end(>=baseN)]
                if d0 in tether_of_buoy:
                    i_buoy, i_dist = d0, d1
                    gidx = tether_of_buoy[i_buoy]
                elif d1 in tether_of_buoy:
                    i_buoy, i_dist = d1, d0
                    gidx = tether_of_buoy[i_buoy]
                else:
                    continue

            # BUOY_LINE: i_dist와 연결된 라인 중 한쪽이 i_dist, 다른쪽이 i_buoy_top
            i_buoy_top = None
            for bl in by_group.get("BUOY_LINE", []):
                n0, n1 = map(int, bl["nodes"])
                if n0 == i_dist:
                    i_buoy_top = n1; break
                if n1 == i_dist:
                    i_buoy_top = n0; break
            if i_buoy_top is None:
                continue

            # Anchor 우선: A/B
            a_seg = None
            for a_candidate in by_group.get("ANCHOR_LINE_A", []):
                n0, n1 = map(int, a_candidate["nodes"])
                if n0 == i_dist:
                    a_seg = (n0, n1); break
                if n1 == i_dist:   # 역순 저장 방어
                    a_seg = (n1, n0); break

            i_anchor_sp, i_anchor_bot = None, None
            if a_seg is not None:
                i_anchor_sp = a_seg[1]
                # B 찾기: 시작이 split이거나 끝이 split
                for b_candidate in by_group.get("ANCHOR_LINE_B", []):
                    n0, n1 = map(int, b_candidate["nodes"])
                    if n0 == i_anchor_sp:
                        i_anchor_bot = n1; break
                    if n1 == i_anchor_sp:
                        i_anchor_bot = n0; break

            # 폴백: 단일 ANCHOR_LINE
            single_anchor = None
            if i_anchor_sp is None or i_anchor_bot is None:
                for al in by_group.get("ANCHOR_LINE", []):
                    n0, n1 = map(int, al["nodes"])
                    if n0 == i_dist:
                        single_anchor = n1; break
                    if n1 == i_dist:
                        single_anchor = n0; break

            # 좌표 계산
            p_base     = coord_of(gidx)
            p_buoy     = coord_of(i_buoy)
            p_dist     = coord_of(i_dist)
            p_buoytop  = coord_of(i_buoy_top)

            if i_anchor_sp is not None and i_anchor_bot is not None:
                p_as = coord_of(i_anchor_sp)
                p_ab = coord_of(i_anchor_bot)
                la = seglen3D(p_dist, p_as)
                lb = seglen3D(p_as,  p_ab)
                a_int, b_int = best_ratio(la, lb, max_den=12)
                p_anchor = p_ab
            elif single_anchor is not None:
                p_as = None
                p_anchor = coord_of(single_anchor)
                a_int, b_int = 1, 1
            else:
                # 앵커 못 찾으면 세트 불완전 → 스킵
                continue

            # 파라미터 역산
            buoy_up     = p_buoy[2]    - p_base[2]
            dist_len    = seglen3D(p_base, p_dist)         # 3D 길이(=수평)
            buoyline_up = p_buoytop[2] - p_dist[2]
            anchor_down = p_dist[2]    - p_anchor[2]

            # 방위각(수평)
            vx, vy = (p_dist[0] - p_base[0], p_dist[1] - p_base[1])
            azimuth_deg = math.degrees(math.atan2(vy, vx)) if (abs(vx) + abs(vy) > 1e-12) else 0.0

            # 경사 각도 역산(수직 기준)
            r_h = seglen3D(p_dist, (p_anchor[0], p_anchor[1], p_dist[2]))
            dz  = max(1e-12, p_dist[2] - p_anchor[2])
            anchor_incl_deg = math.degrees(math.atan(r_h / dz))
            # ✅ Buoy2 사용 여부(새 필드): 길이로 추정
            use_buoy2 = bool(dist_len > 1e-9 and buoyline_up > 1e-9)

            it = {
                "gidx": int(gidx),
                "buoy_up": float(buoy_up),
                "dist_len": float(dist_len),
                "buoyline_up": float(buoyline_up),
                "anchor_down": float(anchor_down),
                "azimuth_deg": float(azimuth_deg),
                "anchor_a": int(a_int),
                "anchor_b": int(b_int),
                "anchor_incl_deg": float(anchor_incl_deg),
                "use_buoy2": use_buoy2,
            }
            restored.append(it)

        # 2) Buoy2 없는 형식: DISTANCE_ROPE가 없으면 폴백 시도
        if not restored and not by_group.get("DISTANCE_ROPE"):
            for ala in by_group.get("ANCHOR_LINE_A", []):
                na0, na1 = map(int, ala["nodes"])
                # base gidx(<baseN) ↔ split(>=baseN) 구조만 수용
                if   na0 < baseN and na1 >= baseN: gidx, i_anchor_sp = na0, na1
                elif na1 < baseN and na0 >= baseN: gidx, i_anchor_sp = na1, na0
                else: continue

                # 대응하는 B 찾기 (split 공유)
                i_anchor_bot = None
                for alb in by_group.get("ANCHOR_LINE_B", []):
                    nb0, nb1 = map(int, alb["nodes"])
                    if   nb0 == i_anchor_sp and nb1 >= baseN: i_anchor_bot = nb1; break
                    elif nb1 == i_anchor_sp and nb0 >= baseN: i_anchor_bot = nb0; break
                if i_anchor_bot is None:
                    continue

                # base에 매단 부이(1) 찾기
                i_buoy = None
                for pl in by_group.get("BUOY_TETHER", []):
                    n0, n1 = map(int, pl["nodes"])
                    if   n0 == gidx and n1 >= baseN: i_buoy = n1; break
                    elif n1 == gidx and n0 >= baseN: i_buoy = n0; break
                if i_buoy is None:
                    continue

                # 좌표
                p_base = coord_of(gidx)
                p_buoy = coord_of(i_buoy)
                p_as   = coord_of(i_anchor_sp)
                p_ab   = coord_of(i_anchor_bot)

                # 파라미터 역산 (Buoy2 없음)
                buoy_up     = p_buoy[2] - p_base[2]
                dist_len    = 0.0
                buoyline_up = 0.0
                anchor_down = p_base[2] - p_ab[2]

                # 방위각: base → anchor_bottom
                vx, vy = (p_ab[0] - p_base[0], p_ab[1] - p_base[1])
                azimuth_deg = math.degrees(math.atan2(vy, vx)) if (abs(vx)+abs(vy))>1e-12 else 0.0

                # 분할 비율로 A:B 근사
                la = seglen3D(p_base, p_as)
                lb = seglen3D(p_as,  p_ab)
                a_int, b_int = best_ratio(la, lb, max_den=12)

                # 경사 각도(수직 기준)
                r_h = ((p_ab[0]-p_base[0])**2 + (p_ab[1]-p_base[1])**2) ** 0.5
                dz  = max(1e-12, p_base[2] - p_ab[2])
                anchor_incl_deg = math.degrees(math.atan(r_h / dz))

                restored.append({
                    "gidx": int(gidx),
                    "buoy_up": float(buoy_up),
                    "dist_len": float(dist_len),
                    "buoyline_up": float(buoyline_up),
                    "anchor_down": float(anchor_down),
                    "azimuth_deg": float(azimuth_deg),
                    "anchor_a": int(a_int),
                    "anchor_b": int(b_int),
                    "anchor_incl_deg": float(anchor_incl_deg),
                    "use_buoy2": False,      # ← Buoy2 없음 명확히 지정
                })

        if not restored:
            return

        # ----- UI 반영 -----
        self.items = restored[:]  # 덮어쓰기
        self.sel_list.clear()

        # 헬퍼가 있으면 쓰고, 없으면 간단 포맷
        if hasattr(self, "_item_label"):
            for it in self.items:
                self.sel_list.addItem(self._item_label(it))
        else:
            for it in self.items:
                label = (f"{it['gidx']} | bu:{it['buoy_up']:.2f}  dist:{it['dist_len']:.2f}  "
                        f"top:{it['buoyline_up']:.2f}  anc:{it['anchor_down']:.2f}  "
                        f"az:{it['azimuth_deg']:.1f}  A:B={it['anchor_a']}:{it['anchor_b']}")
                self.sel_list.addItem(label)

        self.display_vtk()
        print(f"🔁 Reloaded {len(self.items)} item(s) from MOORING_RIG.json (A/B supported).")


    def _preview_segments_for(self, gidx: int):
        import math
        bx, by, bz = self._get_point_by_gidx(int(gidx))

        # 파라미터
        try:
            buoy_up     = float(self.in_buoy_up.text() or 0.0)
            dist_len    = float(self.in_dist_len.text() or 0.0)
            buoyline_up = float(self.in_buoyline_up.text() or 0.0)
            anchor_down = float(self.in_anchor_down.text() or 0.0)
            anchor_deg = float(self.anchor_incl_deg.text() or 0.0)
            a = max(0, int(self.in_anchor_a.text() or "1"))
            b = max(0, int(self.in_anchor_b.text() or "1"))
            if a == 0 and b == 0: a, b = 1, 1
        except Exception:
            buoy_up = dist_len = buoyline_up = anchor_down = 0.0
            a = b = 1

        az_text = (self.in_azimuth_deg.text() or "").strip()
        if az_text:
            try:
                rad = math.radians(float(az_text))
                dx, dy = math.cos(rad), math.sin(rad)
            except Exception:
                dx, dy = 1.0, 0.0
        else:
            n = (bx*bx + by*by) ** 0.5
            dx, dy = ((bx / n) if n > 0 else 1.0, (by / n) if n > 0 else 0.0)

        # 포인트
        p_base      = (bx, by, bz)
        p_buoy      = (bx, by, bz + buoy_up)

        use_buoy2 = bool(getattr(self, "cb_buoy2", None) and self.cb_buoy2.isChecked())
        # Buoy2를 끄면 DISTANCE_ROPE를 0으로 두고 앵커는 p_base에서 바로 시작
        p_dist_end  = (bx + dx * (dist_len if use_buoy2 else 0.0), by + dy * (dist_len if use_buoy2 else 0.0), bz)
        p_buoy_top  = (p_dist_end[0], p_dist_end[1], p_dist_end[2] + (buoyline_up if use_buoy2 else 0.0))

        alpha = math.radians(float(anchor_deg))
        dz = max(0.0, float(anchor_down))            # 수직 낙차 (m)
        r  = dz * math.tan(alpha)                    # 수평 투영거리 (m)
        start_p = p_dist_end if use_buoy2 else p_base
        p_anchor = (start_p[0] + dx * r, start_p[1] + dy * r, start_p[2] - dz)

        # A:B 분할점
        fa = a / (a + b) if (a + b) > 0 else 0.5
        p_split = (
            start_p[0] + fa * (p_anchor[0] - start_p[0]),
            start_p[1] + fa * (p_anchor[1] - start_p[1]),
            start_p[2] + fa * (p_anchor[2] - start_p[2]),
        )

        # 반환 리스트(항상 5개 순서 유지; Buoy2 미사용 시 앞의 2개는 길이 0)
        return [
            (p_base,     p_buoy),      # BUOY_TETHER
            (p_base,     p_dist_end),  # DISTANCE_ROPE (use_buoy2=False면 길이 0)
            (p_dist_end, p_buoy_top),  # BUOY_LINE     (use_buoy2=False면 길이 0)
            (start_p,    p_split),     # ANCHOR_LINE_A
            (p_split,    p_anchor),    # ANCHOR_LINE_B
        ]

    def _collect_item(self, gidx:int):
        buoy_up     = float(self.in_buoy_up.text())
        dist_len    = float(self.in_dist_len.text())
        buoyline_up = float(self.in_buoyline_up.text())
        anchor_down = float(self.in_anchor_down.text())
        az_text = self.in_azimuth_deg.text().strip()
        azimuth_deg = float(az_text) if az_text else None
        try:
            a = max(0, int(self.in_anchor_a.text()))
            b = max(0, int(self.in_anchor_b.text()))
            if a == 0 and b == 0: a, b = 1, 1
        except:
            a, b = 1, 1
        return {
            "gidx": gidx, "buoy_up": buoy_up, "dist_len": dist_len,
            "buoyline_up": buoyline_up, "anchor_down": anchor_down,
            "azimuth_deg": azimuth_deg,
            "anchor_a": a, "anchor_b": b,
            "anchor_incl_deg": self._get_anchor_incl_deg(),\
            "use_buoy2": bool(getattr(self, "cb_buoy2", None) and self.cb_buoy2.isChecked())
        }

    def _preview_segments_for_item(self, it):
        import math
        bx, by, bz = self._get_point_by_gidx(int(it["gidx"]))
        if it.get("azimuth_deg") is not None:
            rad = math.radians(it["azimuth_deg"])
            dx, dy = math.cos(rad), math.sin(rad)
        else:
            n = (bx*bx + by*by) ** 0.5
            dx, dy = ((bx / n) if n > 0 else 1.0, (by / n) if n > 0 else 0.0)

        p_base     = (bx, by, bz)
        p_buoy     = (bx, by, bz + it["buoy_up"])
        use_buoy2 = bool(it.get("use_buoy2", True))
        p_dist_end = (bx + dx * (it["dist_len"] if use_buoy2 else 0.0), by + dy * (it["dist_len"] if use_buoy2 else 0.0), bz)
        p_buoytop  = (p_dist_end[0], p_dist_end[1], p_dist_end[2] + (it["buoyline_up"] if use_buoy2 else 0.0))

        # 앵커 경사: 아이템에 저장된 각도(없으면 현재 UI값) 적용
        alpha = math.radians(float(it.get("anchor_incl_deg", self._get_anchor_incl_deg())))
        dz = max(0.0, float(it["anchor_down"]))
        r = dz * math.tan(alpha)
        start_p = p_dist_end if use_buoy2 else p_base
        p_anchor = (start_p[0] + dx * r, start_p[1] + dy * r, start_p[2] - dz)

        a = int(it.get("anchor_a", 1)); b = int(it.get("anchor_b", 1))
        fa = a / (a + b) if (a + b) > 0 else 0.5
        p_split = (
            start_p[0] + fa * (p_anchor[0] - start_p[0]),
            start_p[1] + fa * (p_anchor[1] - start_p[1]),
            start_p[2] + fa * (p_anchor[2] - start_p[2]),
        )

        return [
            (p_base,     p_buoy),
            (p_base,     p_dist_end),
            (p_dist_end, p_buoytop),
            (start_p,    p_split),   # A
            (p_split,    p_anchor),  # B
         ]

    def build_preview_data(self):
        """
        미리보기는 Mooring frame 노드만 표시.
        display_points: 로컬 인덱스 → 좌표
        display_map: 글로벌 gidx → 로컬 인덱스
        """
        points = [ (float(p[0]), float(p[1]), float(p[2])) for p in self.moor_pts ]
        ids    = [ self.off_moor + i for i in range(len(points)) ]  # ★ 전역 gidx

        # 2) 정규화 + 맵 구성 (문자/정수 혼선 방지용 int(str(...)))
        self.display_points = points
        self.display_ids    = [ int(str(g)) for g in ids ]
        self.display_map    = { self.display_ids[i]: i for i in range(len(self.display_ids)) }

    def display_vtk(self):
        self.renderer.RemoveAllViewProps()

        # =========================
        # 1) Mooring frame 점 (기존 로직)
        # =========================
        pts = vtk.vtkPoints()
        for p in self.display_points:
            pts.InsertNextPoint(p)

        vpoly = vtk.vtkPolyData()
        vpoly.SetPoints(pts)

        verts = vtk.vtkCellArray()
        for i in range(pts.GetNumberOfPoints()):
            verts.InsertNextCell(1)
            verts.InsertCellPoint(i)
        vpoly.SetVerts(verts)

        # 현재/선택 노드 파랑, 나머지 주황
        cur_g = self._current_gidx()

        sel_local = set()
        if cur_g is not None and cur_g in self.display_map:
            sel_local.add(self.display_map[cur_g])
        for gidx in self.sel_nodes:
            if gidx in self.display_map:
                sel_local.add(self.display_map[gidx])

        pcolors = vtk.vtkUnsignedCharArray()
        pcolors.SetNumberOfComponents(3)
        pcolors.SetName("Colors")
        for i in range(len(self.display_points)):
            if i in sel_local:
                pcolors.InsertNextTuple3(51, 128, 255)   # 파랑(선택)
            else:
                pcolors.InsertNextTuple3(230, 100, 70)   # 주황(기본)

        vpoly.GetPointData().SetScalars(pcolors)

        pmap = vtk.vtkPolyDataMapper()
        pmap.SetInputData(vpoly)
        pmap.SetScalarModeToUsePointData()
        pmap.SetColorModeToDirectScalars()
        pmap.ScalarVisibilityOn()

        pactor = vtk.vtkActor()
        pactor.SetMapper(pmap)
        pactor.GetProperty().SetPointSize(7)
        self.renderer.AddActor(pactor)

        # =========================
        # 2) Buoy/Anchor 라인 (그룹별 색)
        # =========================
        # 색상 매핑 (원하면 여기만 바꾸면 됨)
        LINE_COLORS = {
            "BUOY_TETHER":   (120, 180, 255),  # 연파랑
            "DISTANCE_ROPE": (200, 200, 200),  # 회색
            "BUOY_LINE":     ( 30, 110, 255),  # 파랑
            "ANCHOR_LINE_A": (  0, 190,   0),  # 초록
            "ANCHOR_LINE_B": (255, 160,   0),  # 주황
        }
        ORDER = ["BUOY_TETHER", "DISTANCE_ROPE", "BUOY_LINE", "ANCHOR_LINE_A", "ANCHOR_LINE_B"]

        line_pts  = vtk.vtkPoints()
        line_cell = vtk.vtkCellArray()
        lcolors   = vtk.vtkUnsignedCharArray()
        lcolors.SetNumberOfComponents(3)
        lcolors.SetName("Colors")

        # 3) 특별 노드(점) 표시: Buoy1/Buoy2=파랑, Anchor=빨강
        node_pts   = vtk.vtkPoints()
        node_verts = vtk.vtkCellArray()
        ncolors    = vtk.vtkUnsignedCharArray()
        ncolors.SetNumberOfComponents(3)
        ncolors.SetName("Colors")

        def add_mark_point(p, rgb):
            i = node_pts.InsertNextPoint(*p)
            node_verts.InsertNextCell(1)
            node_verts.InsertCellPoint(i)
            ncolors.InsertNextTuple3(*rgb)

        # 대상: 누적 items가 있으면 전부, 없으면 현재 콤보 1개
        def add_segments_from(segs):
            # segs: [(p0,p1), ...] 순서 = ORDER
            for idx, (p0, p1) in enumerate(segs):
                tag = ORDER[idx]
                i0 = line_pts.InsertNextPoint(*p0)
                i1 = line_pts.InsertNextPoint(*p1)
                ln = vtk.vtkLine()
                ln.GetPointIds().SetId(0, i0)
                ln.GetPointIds().SetId(1, i1)
                line_cell.InsertNextCell(ln)
                r, g, b = LINE_COLORS.get(tag, (255, 255, 0))  # 미정 태그는 노랑
                lcolors.InsertNextTuple3(r, g, b)

                # 노드 마킹: Buoy1(p_buoy)=세그 0의 끝점 p1, Buoy2(p_buoy_top)=세그 2의 끝점 p1, Anchor(p_anchor)=세그 4의 끝점 p1
                if idx == 0:  # BUOY_TETHER 끝점 = 부이(1)
                    add_mark_point(p1, (30, 110, 255))   # 파랑
                elif idx == 2:  # BUOY_LINE 끝점 = 부이(2)
                    add_mark_point(p1, (30, 110, 255))   # 파랑
                elif idx == 4:  # ANCHOR_LINE_B 끝점 = 앵커
                    add_mark_point(p1, (220, 40, 40))    # 빨강

        if self.items:
            for it in self.items:
                segs = self._preview_segments_for_item(it)
                add_segments_from(segs)
        else:
            targets = []
            if cur_g is not None:
                targets = [int(cur_g)]
            for g in targets:
                segs = self._preview_segments_for(g)
                add_segments_from(segs)

        # 라인 PolyData (셀별 색상)
        if line_pts.GetNumberOfPoints() > 0:
            lpoly = vtk.vtkPolyData()
            lpoly.SetPoints(line_pts)
            lpoly.SetLines(line_cell)
            lpoly.GetCellData().SetScalars(lcolors)

            lmap = vtk.vtkPolyDataMapper()
            lmap.SetInputData(lpoly)
            lmap.SetScalarModeToUseCellData()
            lmap.SetColorModeToDirectScalars()
            lmap.ScalarVisibilityOn()

            lactor = vtk.vtkActor()
            lactor.SetMapper(lmap)
            lactor.GetProperty().SetLineWidth(2.0)
            self.renderer.AddActor(lactor)

        # Buoy/Anchor 마커 (포인트별 색상)
        if node_pts.GetNumberOfPoints() > 0:
            npoly = vtk.vtkPolyData()
            npoly.SetPoints(node_pts)
            npoly.SetVerts(node_verts)
            npoly.GetPointData().SetScalars(ncolors)

            nmap = vtk.vtkPolyDataMapper()
            nmap.SetInputData(npoly)
            nmap.SetScalarModeToUsePointData()
            nmap.SetColorModeToDirectScalars()
            nmap.ScalarVisibilityOn()

            nactor = vtk.vtkActor()
            nactor.SetMapper(nmap)
            nactor.GetProperty().SetPointSize(9)
            self.renderer.AddActor(nactor)

        # 마무리
        self.renderer.ResetCamera()
        self.renderer.ResetCameraClippingRange()
        self.vtkWidget.GetRenderWindow().Render()

    def add_current_node(self):
        """현재 콤보의 노드와 '입력창의 파라미터'를 한 세트로 캡처해 여러 개 누적."""
        gidx = self.node_combo.currentData()
        if gidx is None:
            return
        it = self._collect_item(int(gidx))         # ← 파라미터 스냅샷
        self.items.append(it)                       # ← 같은 노드도 여러 번 허용
        # 리스트 표시: gidx와 파라미터 요약
        label = (f"{it['gidx']} | bu:{it['buoy_up']}  dist:{it['dist_len']}  "
                f"top:{it['buoyline_up']}  anc:{it['anchor_down']}  "
                f"az:{it['azimuth_deg'] if it['azimuth_deg'] is not None else 'radial'}  "
                f"A:B={it['anchor_a']}:{it['anchor_b']}")
        self.sel_list.addItem(label)
        self.display_vtk()

    def del_selected_node(self):
        """선택된 항목을 self.items에서 제거."""
        row = self.sel_list.currentRow()
        if row >= 0:
            self.sel_list.takeItem(row)
            del self.items[row]
            self.display_vtk()

    # ---------- Data ----------
    def load_data(self):
        """
        현재 장면 기준 전역 인덱스 체계를 복원
        baseN = len(net_pts)+len(float_pts)+len(bc_new_pts)+len(moor_pts)
        """
        with open("Fish_Cage/net_points.json", "r", encoding="utf-8") as f:
            net = json.load(f)
        self.net_pts = net.get("points", [])
        self.z_float = float(net.get("z_float", 0.0))
        self.z_weight = float(net.get("z_weight", -10.0))

        try:
            with open("Fish_Cage/float_temp.json", "r", encoding="utf-8") as f:
                flo = json.load(f)
            self.float_pts = flo.get("points", [])
        except Exception:
            self.float_pts = []

        try:
            with open("Fish_Cage/bottom_collar_temp.json", "r", encoding="utf-8") as f:
                bc = json.load(f)
            # bottom_collar는 points_new만 전역 추가 포인트로 가짐
            self.bc_new_pts = bc.get("points_new", [])
        except Exception:
            self.bc_new_pts = []

        with open("Fish_Cage/mooring_temp.json", "r", encoding="utf-8") as f:
            moor = json.load(f)
        self.moor_pts = moor.get("points", [])

        # offsets
        self.off_net   = 0
        self.off_float = self.off_net + len(self.net_pts)
        self.off_bcnew = self.off_float + len(self.float_pts)
        self.off_moor  = self.off_bcnew + len(self.bc_new_pts)

    def populate_nodes(self):
        self.node_combo.clear()
        for k, p in enumerate(self.moor_pts):
            gi = self.off_moor + k  # 전역 인덱스
            self.node_combo.addItem(f"{gi}: ({p[0]:.2f}, {p[1]:.2f}, {p[2]:.2f})", gi)

    # ---------- Geometry helpers ----------
    @staticmethod
    def _norm2d(vx, vy):
        n = (vx**2 + vy**2) ** 0.5
        return (vx / n, vy / n) if n > 0 else (1.0, 0.0)

    def _centered_radial_dir(self, px, py):
        # 중심에서 바깥으로 (px,py) 방향
        return self._norm2d(px, py)

    def _angle_dir(self, deg):
        rad = math.radians(deg)
        return (math.cos(rad), math.sin(rad))

    def _get_point_by_gidx(self, gidx: int):
        g = int(str(gidx))  # ★ 어떤 타입이 와도 정수화
        if g in self.display_map:
            return self.display_points[self.display_map[g]]
        # fallback: 혹시 gidx가 로컬 인덱스일 수도 있을 때
        if 0 <= g < len(self.display_points):
            return self.display_points[g]
        raise KeyError(f"Unknown gidx: {gidx}")

    # ---------- Save ----------
    def save_rig(self):
        """self.items에 누적된 '노드+파라미터 세트'를 모두 저장.
        비어 있으면 현재 콤보와 현재 입력값으로 1세트 저장."""
        try:
            baseN = len(self.net_pts) + len(self.float_pts) + len(self.bc_new_pts) + len(self.moor_pts)
            next_idx = baseN

            all_extra = []
            all_edges = []
            all_polys = []

            # 저장 대상: items(여러 세트)만 사용 (비어 있으면 '완전 삭제')
            targets = self.items[:]
            if not targets:
                ensure_links_dir()
                out_path = os.path.join(LINKS_DIR, "MOORING_RIG.json")
                if os.path.exists(out_path):
                    try:
                        os.remove(out_path)
                        print(f"🗑️ removed {out_path}")
                    except Exception as e:
                        print(f"⚠️ remove failed {out_path}: {e}")
                # 메인 갱신 후 종료
                if self.parent() and hasattr(self.parent(), "rebuild_main_scene"):
                    self.parent().rebuild_main_scene()
                self.close()
                return

            for it in targets:
                bx, by, bz = self._get_point_by_gidx(int(it["gidx"]))

                # 방향 벡터
                if it.get("azimuth_deg") is not None:
                    rad = math.radians(it["azimuth_deg"])
                    dx, dy = math.cos(rad), math.sin(rad)
                else:
                    n = (bx*bx + by*by) ** 0.5
                    dx, dy = ((bx / n) if n > 0 else 1.0, (by / n) if n > 0 else 0.0)

                # 포인트 생성 (Buoy2 on/off 지원)
                use_buoy2 = bool(it.get("use_buoy2", True))
                buoy_pt   = [bx, by, bz + it["buoy_up"]]

                if use_buoy2:
                    dist_end = [bx + dx * it["dist_len"], by + dy * it["dist_len"], bz]
                    buoy_top = [dist_end[0], dist_end[1], dist_end[2] + it["buoyline_up"]]
                    alpha = math.radians(float(it.get("anchor_incl_deg", self._get_anchor_incl_deg())))
                    dz = max(0.0, float(it["anchor_down"]))
                    r  = dz * math.tan(alpha)
                    anchor_bottom = [dist_end[0] + dx * r, dist_end[1] + dy * r, dist_end[2] - dz]

                    # --- A:B 분할점 계산 ---
                    a = int(it.get("anchor_a", 1)); b = int(it.get("anchor_b", 1))
                    fa = a / (a + b) if (a + b) > 0 else 0.5
                    anchor_split = [
                        dist_end[0] + fa * (anchor_bottom[0] - dist_end[0]),
                        dist_end[1] + fa * (anchor_bottom[1] - dist_end[1]),
                        dist_end[2] + fa * (anchor_bottom[2] - dist_end[2]),
                    ]

                    # --- 인덱스 ---
                    i_buoy       = next_idx; next_idx += 1
                    i_dist_end   = next_idx; next_idx += 1
                    i_buoy_top   = next_idx; next_idx += 1
                    i_anchor_sp  = next_idx; next_idx += 1
                    i_anchor_bot = next_idx; next_idx += 1

                    all_extra.extend([buoy_pt, dist_end, buoy_top, anchor_split, anchor_bottom])

                    # --- 그룹 ---
                    groups = [
                        ("BUOY_TETHER",   [it["gidx"],  i_buoy]),
                        ("DISTANCE_ROPE", [it["gidx"],  i_dist_end]),
                        ("BUOY_LINE",     [i_dist_end,  i_buoy_top]),
                        ("ANCHOR_LINE_A", [i_dist_end,  i_anchor_sp]),  # A
                        ("ANCHOR_LINE_B", [i_anchor_sp, i_anchor_bot]), # B
                    ]
                else:
                    # Buoy2 미사용: DISTANCE_ROPE / BUOY_LINE 생략, 앵커는 base에서 시작
                    alpha = math.radians(float(it.get("anchor_incl_deg", self._get_anchor_incl_deg())))
                    dz = max(0.0, float(it["anchor_down"]))
                    r  = dz * math.tan(alpha)
                    anchor_bottom = [bx + dx * r, by + dy * r, bz - dz]

                    # --- A:B 분할점 계산 ---
                    a = int(it.get("anchor_a", 1)); b = int(it.get("anchor_b", 1))
                    fa = a / (a + b) if (a + b) > 0 else 0.5
                    anchor_split = [
                        bx + fa * (anchor_bottom[0] - bx),
                        by + fa * (anchor_bottom[1] - by),
                        bz + fa * (anchor_bottom[2] - bz),
                    ]

                    # --- 인덱스 ---
                    i_buoy       = next_idx; next_idx += 1
                    i_anchor_sp  = next_idx; next_idx += 1
                    i_anchor_bot = next_idx; next_idx += 1

                    all_extra.extend([buoy_pt, anchor_split, anchor_bottom])

                    # --- 그룹 ---
                    groups = [
                        ("BUOY_TETHER",   [it["gidx"],  i_buoy]),
                        ("ANCHOR_LINE_A", [it["gidx"],  i_anchor_sp]),  # A (base→split)
                        ("ANCHOR_LINE_B", [i_anchor_sp, i_anchor_bot]), # B
                    ]

                for g, nodes in groups:
                    nodes = [int(n) for n in nodes]
                    all_polys.append({"group": g, "nodes": nodes})
                    for a, b in zip(nodes, nodes[1:]):
                        all_edges.append((a, b))

            # 저장
            ensure_links_dir()
            out_path = os.path.join(LINKS_DIR, "MOORING_RIG.json")
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump({
                    "group": "MOORING_RIG",
                    "label": "Mooring rig (buoy, distance, buoy-line, anchor)",
                    "baseN": baseN,
                    "extra_points": all_extra,
                    "edges": all_edges,
                    "polylines": all_polys
                }, f, indent=2, ensure_ascii=False)
            print(f"✅ saved {out_path} (+{len(all_extra)} pts, {len(all_edges)} edges)")

            # 메인 씬 갱신 & 닫기
            if self.parent() and hasattr(self.parent(), "rebuild_main_scene"):
                self.parent().rebuild_main_scene()
            self.close()

        except Exception as e:
            QMessageBox.critical(self, "Save failed", str(e))

class UniversalLinkWindow(QDialog):
    def __init__(self, parent=None, preset=None):
        super().__init__(parent)
        self.preset = preset or {}
        self.setWindowTitle("Link Workbench (Universal)")
        self.resize(960, 680)

        # ---------- UI ----------
        layout = QHBoxLayout(self)

        # VTK view (left)
        self.vtkWidget = QVTKRenderWindowInteractor(self)
        self.renderer = vtk.vtkRenderer()
        self.vtkWidget.GetRenderWindow().AddRenderer(self.renderer)
        self.renderer.SetBackground(0.18, 0.18, 0.18)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.iren.Initialize()
        self.iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        layout.addWidget(self.vtkWidget, stretch=2)

        # Controls (right)
        right = QVBoxLayout()
        right.addWidget(QLabel("Node 1"))
        self.src_section = QComboBox(); right.addWidget(self.src_section)
        self.src_node = QComboBox(); right.addWidget(self.src_node)

        right.addWidget(QLabel("Node 2"))
        self.tgt_section = QComboBox(); right.addWidget(self.tgt_section)
        self.tgt_node = QComboBox(); right.addWidget(self.tgt_node)

        right.addWidget(QLabel("Edge group"))
        self.group_selector = QComboBox()
        self.group_selector.addItem("Floater(braket)", "FLOATER_BRACKET")
        self.group_selector.addItem("Side rope(main)", "SIDE_ROPE_MAIN")
        self.group_selector.addItem("Side rope(sub)",  "SIDE_ROPE_SUB")
        self.group_selector.addItem("Bridle Line", "BRIDLE_LINE")
        right.addWidget(self.group_selector)

        right.addWidget(QLabel("Additional nodes (Side rope main only)"))
        self.addn_input = QLineEdit("0")
        right.addWidget(self.addn_input)

        self.group_selector.currentIndexChanged.connect(self.on_group_changed)

        self.btn_add = QPushButton("Line add"); right.addWidget(self.btn_add)
        self.link_list = QListWidget(); right.addWidget(self.link_list)
        self.btn_del = QPushButton("Line delete"); right.addWidget(self.btn_del)
        self.btn_save = QPushButton("Save"); right.addWidget(self.btn_save)

        layout.addLayout(right, stretch=1)

        # Data
        self.links = []  # list[(src_g, tgt_g, segs, group)]
        self.display_points = []  # preview points
        self.display_map = {}     # global_idx -> local preview idx

        # Signals
        self.src_section.currentIndexChanged.connect(self.reload_nodes)
        self.tgt_section.currentIndexChanged.connect(self.reload_nodes)
        self.btn_add.clicked.connect(self.add_link)
        self.btn_del.clicked.connect(self.del_link)
        self.btn_save.clicked.connect(self.save_links)
        self.src_node.currentIndexChanged.connect(self.display_vtk)
        self.tgt_node.currentIndexChanged.connect(self.display_vtk)

         # --- 데이터 로드 & 섹션 초기화 ---
        self.load_data()
        self.init_sections()   # 여기서 self.src_section/self.tgt_section 에 addItem(key, key) 형태로 들어가야 함

        self._apply_boot_flag()    # 🔴 여기 한 줄 추가
        self.prefill_from_group(self.group_selector.currentData())  # 🔹 현재 그룹 자동 로드
        self.reload_nodes()
        self.display_vtk()
        self.on_group_changed()  # ← 초기 상태 반영
        self._warned_invalid_indices = set()  # 잘못된 인덱스 팝업 중복 방지


    def on_group_changed(self):
        gkey = self.group_selector.currentData()
        # 추가 노드 입력: SIDE_ROPE_MAIN만 사용
        if hasattr(self, "addn_input"):
            self.addn_input.setEnabled(gkey == "SIDE_ROPE_MAIN")

        # BRIDLE_LINE이면 기본 섹션 자동 세팅: out-floater → mooring-frame
        if gkey == "BRIDLE_LINE":
            s = self.src_section.findData("out-floater")
            t = self.tgt_section.findData("mooring-frame")
            if s >= 0: self.src_section.setCurrentIndex(s)
            if t >= 0: self.tgt_section.setCurrentIndex(t)
            # 섹션 바뀌었으니 즉시 노드 리스트 갱신
            self.reload_nodes()
        # ✅ 선택한 그룹의 기존 링크만 로드(세션 간 섞임 방지)
        self.prefill_from_group(gkey)

    # ---------- Data/loading ----------
    def load_data(self):
        with open("Fish_Cage/net_points.json") as f:
            net = json.load(f)

        self.net_pts = net.get("points", [])
        self.z_float = float(net.get("z_float", 0.0))
        self.z_weight = float(net.get("z_weight", -10.0))

        # groups는 더 이상 필수 아님 (없어도 동작), 있으면 호환
        raw_groups = net.get("groups", {}) or {}
        self.net_groups = {
            k: [int(x) for x in (raw_groups.get(k) or [])]
            for k in ("top_ring", "bottom_ring", "bottom_plate")
        }

        # side/bottom JSON: 좌표 or 인덱스 모두 지원
        def _try_load_indices_from_side_bottom(path_key):
            try:
                with open(path_key, "r", encoding="utf-8") as f:
                    d = json.load(f)
                if "points" in d and isinstance(d["points"], list):
                    return self._match_points_to_indices(d["points"])
                idxs = d.get("indices", [])
                return [int(i) for i in (idxs or [])]
            except Exception:
                return []

        self.side_indices_file = _try_load_indices_from_side_bottom("Fish_Cage/side_net.json")
        self.bottom_indices_file = _try_load_indices_from_side_bottom("Fish_Cage/bottom_net.json")

        # Float (optional)
        try:
            with open("Fish_Cage/float_temp.json") as f:
                flo = json.load(f)
            self.float_pts = flo.get("points", [])
        except Exception:
            self.float_pts = []

        # Bottom collar (optional)
        try:
            with open("Fish_Cage/bottom_collar_temp.json") as f:
                bc = json.load(f)
            self.bc_new_pts = bc.get("points_new", [])
            self.bc_ring_idxs = [int(i) for i in bc.get("ring_indices", [])]
        except Exception:
            self.bc_new_pts = []
            self.bc_ring_idxs = []

        # Mooring frame (optional)
        try:
            with open("Fish_Cage/mooring_temp.json", "r", encoding="utf-8") as f:
                moor = json.load(f)
            self.moor_pts = moor.get("points", [])
        except Exception:
            self.moor_pts = []

        # Global offsets (must match main rebuild order)
        self.off_net = 0
        self.off_float = len(self.net_pts)
        self.off_bcnew = self.off_float + len(self.float_pts)
        self.off_moor  = self.off_bcnew + len(self.bc_new_pts)   # ✅ 추가

        # Previous baked universal (optional)
        try:
            with open("Fish_Cage/link_universal_temp.json") as f:
                self.prev_universal = json.load(f)
        except Exception:
            self.prev_universal = {}

        self.link_extra_points = {}  # {global_idx: [x,y,z]}
        try:
            ensure_links_dir()
            for path in sorted(glob.glob(os.path.join(LINKS_DIR, "*.json"))):
                with open(path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                baseN = int(data.get("baseN", 0))
                extra = data.get("extra_points", []) or []
                for k, pt in enumerate(extra):
                    self.link_extra_points[baseN + k] = pt
        except Exception as e:
            print(f"ℹ️ link-extra load skipped: {e}")

    # ---------- Section builders ----------
    def _in_floater_indices(self):
        return [i for i, p in enumerate(self.net_pts) if abs(p[2] - self.z_float) < 1e-6]

    def _out_floater_indices(self):
        return list(range(self.off_float, self.off_float + len(self.float_pts)))

    def _bottom_collar_indices(self):
        return list(self.bc_ring_idxs)
    
    def _mooring_frame_indices(self):
        return list(range(self.off_moor, self.off_moor + len(self.moor_pts)))

    def _side_rope_main_indices(self):
        nodes = []

        # 1) 그룹별 파일 우선
        try:
            path = os.path.join(LINKS_DIR, "SIDE_ROPE_MAIN.json")
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            for poly in (data.get("polylines", []) or []):
                for n in (poly.get("nodes", []) or []):
                    try:
                        nodes.append(int(n))
                    except Exception:
                        pass
        except Exception:
            # 2) 폴백: 합본(universal)
            try:
                path = "Fish_Cage/link_universal_temp.json"
                with open(path, "r", encoding="utf-8") as f:
                    uni = json.load(f)
                for poly in (uni.get("polylines", []) or []):
                    if poly.get("group") == "SIDE_ROPE_MAIN":
                        for n in (poly.get("nodes", []) or []):
                            try:
                                nodes.append(int(n))
                            except Exception:
                                pass
            except Exception:
                nodes = []

        # 중복 제거(등장 순서 보존) + 범위 방어
        seen, out = set(), []
        for gidx in nodes:
            if gidx not in seen:
                out.append(int(gidx)); seen.add(int(gidx))
        return out
    
    def _bottom_net_indices(self):
        if getattr(self, "bottom_indices_file", []):
            n = len(self.net_pts)
            return [i for i in self.bottom_indices_file if 0 <= i < n]

        tol = 1e-6
        idxs = [gi for gi, p in enumerate(self.net_pts) if p[2] <= self.z_weight + tol]
        idxs.sort()
        return idxs

    def _side_net_indices(self):
        if getattr(self, "side_indices_file", []):
            n = len(self.net_pts)
            return [i for i in self.side_indices_file if 0 <= i < n]

        tol = 1e-6
        n = len(self.net_pts)

        top_ring = set((getattr(self, "net_groups", {}) or {}).get("top_ring", []))
        if not top_ring:
            top_ring = {i for i, p in enumerate(self.net_pts) if abs(p[2] - self.z_float) < tol}

        bottom_plate = set((getattr(self, "net_groups", {}) or {}).get("bottom_plate", []))
        if not bottom_plate:
            bottom_plate = {i for i, p in enumerate(self.net_pts) if p[2] <= self.z_weight + tol}

        # 범위 방어
        top_ring = {i for i in top_ring if 0 <= i < n}
        bottom_plate = {i for i in bottom_plate if 0 <= i < n}

        return [i for i in range(n) if i not in top_ring and i not in bottom_plate]


    def build_sections(self):
        return {
            "in-floater": self._in_floater_indices(),
            "out-floater": self._out_floater_indices(),
            "bottom-collar": self._bottom_collar_indices(),
            "side-net": self._side_net_indices(),
            "bottom-net": self._bottom_net_indices(),
            "mooring-frame": self._mooring_frame_indices(),   # ✅ 추가
            "side-rope-main": self._side_rope_main_indices(),
            
        }

    def init_sections(self):
        self.sections = self.build_sections()
        self.src_section.clear(); self.tgt_section.clear()
        for key in ["in-floater", "out-floater", "bottom-collar", "side-net", "bottom-net", "mooring-frame", "side-rope-main"]:
            self.src_section.addItem(key, key)
            self.tgt_section.addItem(key, key)
        # sensible default: out-floater -> bottom-collar
        self.src_section.setCurrentIndex(0)
        self.tgt_section.setCurrentIndex(1)

    def _apply_boot_flag(self):
        """글로벌 BOOT_FLAG에 따라 기본 세팅 강제"""
        global BOOT_FLAG
        if BOOT_FLAG == 0:
            # 0: Bracket (IN→OUT)
            self.src_section.setCurrentIndex(0)
            self.tgt_section.setCurrentIndex(1)
            gi = self.group_selector.findData("FLOATER_BRACKET")
            if gi < 0: gi = self.group_selector.findText("Floater(braket)")
            if gi >= 0: self.group_selector.setCurrentIndex(gi)
            # 예: 세그먼트 1
            # self.seg_input.setText("1")

        elif BOOT_FLAG == 1:
            # 1: Side rope(main) (OUT→BC) — 당신 로직대로 1,2 인덱스 사용
            self.src_section.setCurrentIndex(1)
            self.tgt_section.setCurrentIndex(2)
            # ✅ 그룹: Side rope(main)
            gi = self.group_selector.findData("SIDE_ROPE_MAIN")
            if gi < 0: gi = self.group_selector.findText("Side rope(main)")
            if gi >= 0: self.group_selector.setCurrentIndex(gi)

        elif BOOT_FLAG == 2:
            # 2: Side rope(sub) (OUT→BC)
            self.src_section.setCurrentIndex(4)   # bottom-net
            self.tgt_section.setCurrentIndex(2)   # bottom-collar
            gi = self.group_selector.findData("SIDE_ROPE_SUB")
            if gi < 0: gi = self.group_selector.findText("Side rope(sub)")
            if gi >= 0: self.group_selector.setCurrentIndex(gi)

        elif BOOT_FLAG == 3:
            # 3: Bridle Line (OUT → Mooring)
            gi = self.group_selector.findData("BRIDLE_LINE")
            if gi < 0: gi = self.group_selector.findText("Bridle Line")
            if gi >= 0: self.group_selector.setCurrentIndex(gi)
            s = self.src_section.findData("out-floater")
            t = self.tgt_section.findData("mooring-frame")
            if s >= 0: self.src_section.setCurrentIndex(s)
            if t >= 0: self.tgt_section.setCurrentIndex(t)
        BOOT_FLAG = -1




    # ---------- UI refresh ----------
    def reload_nodes(self):
        self.sections = self.build_sections()
        s_key = self.src_section.currentData()
        t_key = self.tgt_section.currentData()
        s_list = self.sections.get(s_key, [])
        t_list = self.sections.get(t_key, [])

        self.src_node.clear()
        for gi in s_list:
            x, y, z = self.get_point(gi)
            self.src_node.addItem(f"S{gi}: ({x:.2f},{y:.2f},{z:.2f})", gi)

        self.tgt_node.clear()
        for gi in t_list:
            x, y, z = self.get_point(gi)
            self.tgt_node.addItem(f"T{gi}: ({x:.2f},{y:.2f},{z:.2f})", gi)

        # Preview points: unique union of both lists
        uniq, seen = [], set()
        for gi in s_list + t_list:
            if gi not in seen:
                uniq.append(gi); seen.add(gi)
        self.display_points = [self.get_point(g) for g in uniq]
        self.display_map = {g: i for i, g in enumerate(uniq)}
        self.display_order = uniq[:]   # ⭐ 로컬 i → 글로벌 gidx 매핑 저장
        self.display_vtk()

    def get_point(self, gidx):
        # 1) 기본 스택: net → float → bottom-collar → mooring
        if gidx < self.off_float:
            return self.net_pts[gidx]
        elif gidx < self.off_bcnew:
            return self.float_pts[gidx - self.off_float]
        elif gidx < self.off_moor:
            return self.bc_new_pts[gidx - self.off_bcnew]
        else:
            idx = gidx - self.off_moor
            if 0 <= idx < len(self.moor_pts):
                return self.moor_pts[idx]

        # 2) 저장된 추가 노드(링크 extra)를 마지막에 시도
        if hasattr(self, "link_extra_points") and gidx in self.link_extra_points:
            return self.link_extra_points[gidx]

        # 3) 안전망(팝업 경고는 그대로 유지)
        try:
            bad = int(gidx)
        except Exception:
            bad = gidx
        if bad not in getattr(self, "_warned_invalid_indices", set()):
            self._warned_invalid_indices.add(bad)
            QMessageBox.warning(
                self, "Invalid node index",
                f"유효하지 않은 노드 인덱스가 참조되었습니다: {bad}\n"
                f"저장 파일(link/universal)과 현재 장면의 인덱스 기준이 어긋났을 수 있습니다."
            )
        return [0.0, 0.0, 0.0]
    
    # ---------- Link ops ----------
    def add_link(self):
        s = self.src_node.currentData()
        e = self.tgt_node.currentData()
        gkey = self.group_selector.currentData()
        glab = self.group_selector.currentText()

        # ✅ SIDE_ROPE_MAIN만 "추가 노드 개수" 적용
        if gkey == "SIDE_ROPE_MAIN":
            try:
                add_n = max(0, int(self.addn_input.text()))  # 추가 노드 개수(0 이상)
            except Exception:
                add_n = 0
            segs = add_n + 1       # 내부 보간 로직은 segs-1 개의 중간점을 만듦
        else:
            add_n = 0
            segs = 1               # 다른 그룹은 중간점 없음

        if s is None or e is None:
            return
        tup = (int(s), int(e), int(segs), gkey)
        if tup not in self.links:
            self.links.append(tup)
            self.link_list.addItem(
                f"{self.src_section.currentData()}:{s} ↔ {self.tgt_section.currentData()}:{e}  (+nodes:{add_n}, grp:{glab})"
            )
            self.display_vtk()

    def del_link(self):
        row = self.link_list.currentRow()
        if row >= 0:
            self.link_list.takeItem(row)
            del self.links[row]
            self.display_vtk()

    # ---------- Save (bake to JSON) ----------
    def save_links(self):
        ensure_links_dir()
        baseN = len(self.net_pts) + len(self.float_pts) + len(self.bc_new_pts) + len(self.moor_pts)

        # 그룹별 묶기
        grouped = {}
        for (s,e,segs,gkey) in self.links:
            grouped.setdefault(gkey, []).append((int(s), int(e), int(segs)))

        labels_map = {self.group_selector.itemData(i): self.group_selector.itemText(i)
                    for i in range(self.group_selector.count())}

        #   ▶ '현재 선택된 그룹'을 비워서 저장한 경우에만 그 그룹 파일 삭제
        #   ▶ 다른 그룹 파일은 건드리지 않음 (보존)
        current_gkey = self.group_selector.currentData()
        if not grouped.get(current_gkey):  # 현재 그룹 항목이 하나도 없음(=비움)
            ensure_links_dir()
            path = os.path.join(LINKS_DIR, f"{current_gkey}.json")
            if os.path.exists(path):
                os.remove(path)

        all_polylines = []  # 🔹 universal 호환 저장용

        for gkey, items in grouped.items():
            extra_pts, edges, polylines = [], [], []
            next_idx = baseN

            def P(gidx):
                n_net   = len(self.net_pts)
                n_float = len(self.float_pts)
                n_bc    = len(self.bc_new_pts)
                n_moor  = len(self.moor_pts)

                if gidx < n_net:
                    return self.net_pts[gidx]
                elif gidx < n_net + n_float:
                    return self.float_pts[gidx - n_net]
                elif gidx < n_net + n_float + n_bc:
                    off = n_net + n_float
                    return self.bc_new_pts[gidx - off]
                elif gidx < n_net + n_float + n_bc + n_moor:
                    off = n_net + n_float + n_bc
                    return self.moor_pts[gidx - off]
                else:
                    return None

            for (s,e,segs) in items:
                ps, pe = P(s), P(e)
                prev = s
                local_new = []
                if segs > 1:
                    for k in range(1, segs):
                        t = k / segs
                        x = (1-t)*ps[0] + t*pe[0]
                        y = (1-t)*ps[1] + t*pe[1]
                        z = (1-t)*ps[2] + t*pe[2]
                        extra_pts.append([x,y,z])
                        cur = next_idx; next_idx += 1
                        edges.append((prev, cur)); prev = cur
                        local_new.append(cur)
                    edges.append((prev, e))
                else:
                    edges.append((s, e))

                poly = {"group": gkey, "nodes": [s] + local_new + [e]}
                polylines.append(poly)
                all_polylines.append(poly)  # 🔹 universal 합본에도 추가

            path = os.path.join(LINKS_DIR, f"{gkey}.json")
            with open(path, "w", encoding="utf-8") as f:
                json.dump({
                    "group": gkey,
                    "label": labels_map.get(gkey, gkey),
                    "baseN": baseN,
                    "extra_points": extra_pts,
                    "edges": edges,
                    "polylines": polylines
                }, f, indent=2, ensure_ascii=False)
            print(f"✅ saved {path}  (+{len(extra_pts)} pts, {len(edges)} edges)")

            # 다음 그룹도 같은 baseN에서 누적 생성하도록 유지 (각 그룹은 독립 파일)
            baseN = baseN + len(extra_pts)

        # 🔹 호환성: 예전 로더를 위한 합본도 갱신
        with open("Fish_Cage/link_universal_temp.json", "w", encoding="utf-8") as f:
            json.dump({"polylines": all_polylines}, f, indent=2, ensure_ascii=False)
        print("✅ link_universal_temp.json updated")

        # 메인 갱신
        if self.parent() and hasattr(self.parent(), "rebuild_main_scene"):
            self.parent().rebuild_main_scene()
        self.close()


    # ---------- VTK preview ----------
    def display_vtk(self):
        self.renderer.RemoveAllViewProps()

        # Points
        pts = vtk.vtkPoints()
        for p in self.display_points:
            pts.InsertNextPoint(p)

        vpoly = vtk.vtkPolyData()
        vpoly.SetPoints(pts)

        verts = vtk.vtkCellArray()
        for i in range(pts.GetNumberOfPoints()):
            verts.InsertNextCell(1)
            verts.InsertCellPoint(i)
        vpoly.SetVerts(verts)

        # --- 선택 노드 2개만 파란색으로 칠하기 (POINT DATA SCALARS 강제 사용) ---
        s_g = self.src_node.currentData()
        e_g = self.tgt_node.currentData()
        try:
            s_g = None if s_g is None else int(s_g)
            e_g = None if e_g is None else int(e_g)
        except Exception:
            pass

        colors = vtkUnsignedCharArray()
        colors.SetNumberOfComponents(3)
        colors.SetName("Colors")

        sel_local = set()
        if s_g is not None and s_g in self.display_map:
            sel_local.add(self.display_map[s_g])
        if e_g is not None and e_g in self.display_map:
            sel_local.add(self.display_map[e_g])

        # 2) 로컬 인덱스 기준으로 색 지정
        for i in range(len(self.display_points)):
            if i in sel_local:
                colors.InsertNextTuple3(51, 128, 255)   # 파랑(선택)
            else:
                colors.InsertNextTuple3(230, 100, 70)   # 주황(기본)

        vpoly.GetPointData().SetScalars(colors)

        pmap = vtk.vtkPolyDataMapper()
        pmap.SetInputData(vpoly)
        # ⬇️ 스칼라 컬러를 '무조건' 쓰게 강제
        pmap.SetScalarModeToUsePointData()
        pmap.SetColorModeToDirectScalars()
        pmap.ScalarVisibilityOn()

        pactor = vtk.vtkActor()
        pactor.SetMapper(pmap)
        pactor.GetProperty().SetPointSize(6)

        # Lines: preview as straight lines between currently added pairs
        lines = vtk.vtkCellArray()
        for (s, e, _segs, _g) in self.links:
            if s in self.display_map and e in self.display_map:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, self.display_map[s])
                line.GetPointIds().SetId(1, self.display_map[e])
                lines.InsertNextCell(line)

        lpoly = vtk.vtkPolyData(); lpoly.SetPoints(pts); lpoly.SetLines(lines)
        lmap = vtk.vtkPolyDataMapper(); lmap.SetInputData(lpoly)
        lactor = vtk.vtkActor(); lactor.SetMapper(lmap)
        lactor.GetProperty().SetColor(1.0, 1.0, 0.0)

        self.renderer.AddActor(pactor)
        self.renderer.AddActor(lactor)
        self.renderer.ResetCamera()
        self.vtkWidget.GetRenderWindow().Render()

    def apply_preset_safe(self):
        if not self.preset:
            return

        def _set_by_data_or_text(combo: QComboBox, value: str):
            if not value:
                return
            # 1) userData로 찾기
            i = combo.findData(value)
            if i < 0:
                # 2) 텍스트로도 시도 (UI 라벨과 동일한 문자열인 경우)
                j = combo.findText(value)
                if j >= 0:
                    combo.setCurrentIndex(j)
                    return
            if i >= 0:
                combo.setCurrentIndex(i)

        # 섹션 (Node1 / Node2)
        _set_by_data_or_text(self.src_section, self.preset.get("src"))  # 예: "in-floater"
        _set_by_data_or_text(self.tgt_section, self.preset.get("tgt"))  # 예: "out-floater"

        # 그룹
        _set_by_data_or_text(self.group_selector, self.preset.get("group"))  # "FLOATER_BRACKET"

        # 현재 선택 상태로 리스트/프리뷰 갱신
        self.reload_nodes()

        # 저장본 채우기 (있으면)
        try:
            self.prefill_from_group(self.preset.get("group"))
        except Exception as e:
            print(f"[preset] prefill skipped: {e}")

    def prefill_from_group(self, group_key: str):
        if not group_key:
            return

        polys = []
        # 1) links/{group}.json 우선
        try:
            path = os.path.join(LINKS_DIR, f"{group_key}.json")
            with open(path, "r", encoding="utf-8") as f:
                data = json.load(f)
            polys = data.get("polylines", []) or []
        except Exception:
            polys = []

        # 2) 폴백: link_universal_temp.json
        if not polys:
            try:
                with open("Fish_Cage/link_universal_temp.json", "r", encoding="utf-8") as f:
                    uni = json.load(f)
                polys = [p for p in (uni.get("polylines", []) or []) if p.get("group") == group_key]
            except Exception:
                polys = []

        self.links.clear()
        self.link_list.clear()

        glabel = self.group_selector.currentText()
        loaded = 0
        for poly in polys:
            nodes = poly.get("nodes", [])
            if len(nodes) < 2: 
                continue
            s, e = int(nodes[0]), int(nodes[-1])
            segs = max(1, len(nodes) - 1)
            self.links.append((s, e, segs, group_key))
            self.link_list.addItem(
                f"{self.src_section.currentData()}:{s} ↔ {self.tgt_section.currentData()}:{e} (seg:{segs}, grp:{glabel})"
            )
            loaded += 1

        self.display_vtk()
        print(f"ℹ️ Loaded {loaded} ropes for group={group_key}")
        
    def _match_points_to_indices(self, pts_list, ndigits=6):
        """
        좌표 리스트 → self.net_pts 인덱스로 매칭
        """
        def _key(p):
            return (round(float(p[0]), ndigits), round(float(p[1]), ndigits), round(float(p[2]), ndigits))

        lookup = {}
        for gi, p in enumerate(self.net_pts):
            lookup.setdefault(_key(p), []).append(gi)

        out = []
        for p in pts_list:
            key = _key(p)
            if key in lookup and lookup[key]:
                gi = lookup[key].pop(0)
                out.append(gi)
        return sorted(set(out))
    
class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setWindowTitle("Fish Cage Designer (DaeYeon Cho)")
        self.resize(1200, 800)

        frame = QFrame()
        self.setCentralWidget(frame)
        h_layout = QHBoxLayout()

        self.vtkWidget = QVTKRenderWindowInteractor(frame)

        self.renderer = vtk.vtkRenderer()
        self.vtkWidget.GetRenderWindow().AddRenderer(self.renderer)
        self.renderer.SetBackground(0.1, 0.1, 0.1)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()
        self.iren.Initialize()
        self.iren.SetInteractorStyle(PanOnlyStyle())

        # ── 메뉴: Views → X axis / Y axis / Z axis ──────────────────────
        views = self.menuBar().addMenu("&Views")
        act_x = QAction("X axis", self); act_x.setToolTip("View along +X  (Shift: -X)")
        act_y = QAction("Y axis", self); act_y.setToolTip("View along +Y  (Shift: -Y)")
        act_z = QAction("Z axis", self); act_z.setToolTip("View along +Z  (Shift: -Z)")
        act_x.triggered.connect(lambda: self._view_axis('x'))
        act_y.triggered.connect(lambda: self._view_axis('y'))
        act_z.triggered.connect(lambda: self._view_axis('z'))
        views.addActions([act_x, act_y, act_z])
         # ────────────────────────────────────────────────────────────────
        # ── 메뉴: Data (JSON 삭제/초기화) ───────────────────────────────
        data_menu = self.menuBar().addMenu("&Delete Data")

        act_del_float    = QAction("Delete Floater Data", self)
        act_del_bottom   = QAction("Delete bottom collar", self)
        act_del_moor     = QAction("Delete mooring Frame", self)
        act_del_mrig     = QAction("Delete MOORING Line", self)
        act_del_links    = QAction("Delete links/*.json (all link groups)", self)
        act_del_uni      = QAction("Delete link_universal_temp.json", self)
        act_del_legacy1  = QAction("Delete floating_link_temp.json", self)
        act_del_legacy2  = QAction("Delete side_rope_link_temp.json", self)

        act_reset_all    = QAction("✅ 전체 초기화 (Fish_Cage 내 JSON 전부 삭제)", self)

        data_menu.addActions([
            act_del_float, act_del_bottom,
            act_del_moor, act_del_mrig, act_del_links, act_del_uni,
            act_del_legacy1, act_del_legacy2
        ])
        data_menu.addSeparator()
        data_menu.addAction(act_reset_all)

        act_del_float.triggered.connect(lambda: self._delete_confirm([
            "Fish_Cage/float_temp.json"
        ], "float_temp.json을 삭제할까요?"))

        act_del_bottom.triggered.connect(lambda: self._delete_confirm([
            "Fish_Cage/bottom_collar_temp.json"
        ], "bottom_collar_temp.json을 삭제할까요?"))

        act_del_moor.triggered.connect(lambda: self._delete_confirm([
            "Fish_Cage/mooring_temp.json"
        ], "mooring_temp.json을 삭제할까요?"))

        act_del_mrig.triggered.connect(lambda: self._delete_confirm([
            "Fish_Cage/links/MOORING_RIG.json"
        ], "links/MOORING_RIG.json을 삭제할까요?"))

        act_del_links.triggered.connect(lambda: self._delete_links_dir("Fish_Cage/links"))

        act_del_uni.triggered.connect(lambda: self._delete_confirm([
            "Fish_Cage/link_universal_temp.json"
        ], "link_universal_temp.json을 삭제할까요?"))

        act_del_legacy1.triggered.connect(lambda: self._delete_confirm([
            "Fish_Cage/floating_link_temp.json"
        ], "floating_link_temp.json을 삭제할까요?"))

        act_del_legacy2.triggered.connect(lambda: self._delete_confirm([
            "Fish_Cage/side_rope_link_temp.json"
        ], "side_rope_link_temp.json을 삭제할까요?"))

        act_reset_all.triggered.connect(self._full_reset_fish_cage)
        # ────────────────────────────────────────────────────────────────

        self._install_axes_gizmo(self.iren)

        self.tab_widget = QTabWidget()

        # ▶ 기존 '왼쪽 VTK + 오른쪽 패널'을 QSplitter로 감싼다
        self.splitter = QSplitter(Qt.Horizontal, frame)
        self.splitter.setObjectName("main_splitter")  # 스타일 범위 지정
        self.splitter.addWidget(self.vtkWidget)     # 왼쪽
        self.splitter.addWidget(self.tab_widget)    # 오른쪽
        self.splitter.setStretchFactor(0, 1)        # 왼쪽 가변
        self.splitter.setStretchFactor(1, 0)        # 오른쪽 고정 느낌
        # ▶ 손잡이 폭을 넉넉히 (버튼이 잘리지 않게)
        #   * 테마/OS 기본값(5~9px)은 26x26 버튼이 대부분 안 보입니다.
        # ▶ 손잡이는 충분히 넓히되, '완전 투명'으로 만들어 배경과 동일하게 보이게
        #    (버튼만 보이도록)
        self.splitter.setHandleWidth(10)
        self.splitter.setSizes([900, 360])
        h_layout.addWidget(self.splitter)

        # ▶ 손잡이에 화살표 토글 버튼(>>, <<) 설치
        self._install_right_panel_toggle(self.splitter, handle_index=1, default_right=320)

        self.tab_net = QWidget()
        form_net = QFormLayout()

        self.sides_input = QLineEdit("24")
        self.f_circum_input = QLineEdit("100")
        self.z_float_input = QLineEdit("0")
        self.bn_circum_input = QLineEdit("100")
        self.z_weight_input = QLineEdit("-30")
        self.half_mesh_input = QLineEdit("1.0")
        self.diameter_input = QLineEdit("0.3")
        self.bottomnetangle_input = QLineEdit("15")

        form_net.addRow(QLabel("Sides"), self.sides_input)
        form_net.addRow(QLabel("Net Circumference (m)"), self.f_circum_input)
        form_net.addRow(QLabel("Z Position (m)"), self.z_float_input)

        hline1 = QFrame()
        hline1.setFrameShape(QFrame.HLine)
        hline1.setFrameShadow(QFrame.Sunken)
        form_net.addRow(hline1)

        form_net.addRow(QLabel("Bottom Net Circumference (m)"), self.bn_circum_input)
        form_net.addRow(QLabel("Z Position (m)"), self.z_weight_input)
        
        geom_btn = QPushButton("Submit Geometry")
        geom_btn.clicked.connect(self.submit_geometry)
        form_net.addRow(geom_btn)

        hline1 = QFrame()
        hline1.setFrameShape(QFrame.HLine)
        hline1.setFrameShadow(QFrame.Sunken)
        form_net.addRow(hline1)

        form_net.addRow(QLabel("Half Mesh Size (m)"), self.half_mesh_input)
        form_net.addRow(QLabel("Twine Diameter (m)"), self.diameter_input)
        form_net.addRow(QLabel("BottomNet Angle (°)"), self.bottomnetangle_input)

        net_btn = QPushButton("Submit Net Parameter")
        net_btn.clicked.connect(self.submit_net_parameter)
        form_net.addRow(net_btn)

        hline1 = QFrame()
        hline1.setFrameShape(QFrame.HLine)
        hline1.setFrameShadow(QFrame.Sunken)
        form_net.addRow(hline1)

        form_net.addRow(self.create_common_buttons())

        self.tab_net.setLayout(form_net)
        self.tab_widget.addTab(self.tab_net, "Net")

###########################################################################
#                               Collar 탭                                 #
###########################################################################
        self.tab_collar = QWidget()
        form_collar = QFormLayout()

        title_label = QLabel("Floating Collar")
        title_label.setAlignment(Qt.AlignCenter)
        form_collar.addRow(title_label)  

        self.bracket_num = QLineEdit("0")
        self.float_outside_cir = QLineEdit("105")
        self.outside_sides = QLineEdit("24")


        form_collar.addRow(QLabel("Number of bracket"), self.bracket_num)
        form_collar.addRow(QLabel("Outside Floater Sides"), self.outside_sides)        
        form_collar.addRow(QLabel("circumference (m)"), self.float_outside_cir)

        self.float_generate_btn = QPushButton("Generate Outer Floater")
        self.float_generate_btn.clicked.connect(self.generate_outer_floater_json)
        form_collar.addRow(self.float_generate_btn)

        self.open_floating_link_btn = QPushButton("Open Floating Bracket Link")
        self.open_floating_link_btn.clicked.connect(self.open_floating_link_window)
        form_collar.addRow(self.open_floating_link_btn)

        hline1 = QFrame()
        hline1.setFrameShape(QFrame.HLine)
        hline1.setFrameShadow(QFrame.Sunken)
        form_collar.addRow(hline1)

        title_label = QLabel("Bottom Collar")
        title_label.setAlignment(Qt.AlignCenter)
        form_collar.addRow(title_label) 


        # ✅ 모드 선택 (기본: Sinker array)
        self.rope_mode_group = QButtonGroup(self)
        self.rb_sinkers = QRadioButton("Sinker array (one per rope)")
        self.rb_bottom_collar = QRadioButton("Bottom collar (ring)")
        self.rb_sinkers.setChecked(True)
        self.rope_mode_group.addButton(self.rb_sinkers)
        self.rope_mode_group.addButton(self.rb_bottom_collar)
        form_collar.addRow(self.rb_sinkers)
        form_collar.addRow(self.rb_bottom_collar)

        self.bottomcollar_cir = QLineEdit("")
        self.bottomcollar_z = QLineEdit("")

        # ✅ circumference 라벨을 변수로 잡아둠 (enable/disable 용)
        self.lbl_bottomcollar_cir = QLabel("circumference (m)")
        self.lbl_bottomcollar_z = QLabel("Z position (m)")
        self.btn_gen_bottom_collar = QPushButton("Generate Bottom Collar")
        self.open_side_rope_link_btn = QPushButton("Open Side Rope Link (main)")
        self.open_side_rope_link_sub_btn = QPushButton("Open Side Rope Link (Sub)")

        form_collar.addRow(self.lbl_bottomcollar_cir, self.bottomcollar_cir)
        form_collar.addRow(self.lbl_bottomcollar_z, self.bottomcollar_z)
        self.btn_gen_bottom_collar.clicked.connect(self.generate_bottom_collar)
        form_collar.addRow(self.btn_gen_bottom_collar)

        self.open_side_rope_link_btn.clicked.connect(self.open_side_rope_link_window)
        form_collar.addRow(self.open_side_rope_link_btn)
        self.open_side_rope_link_sub_btn.clicked.connect(self.open_side_rope_link_sub_window)
        form_collar.addRow(self.open_side_rope_link_sub_btn)

        # ✅ 토글 시 핸들러 연결
        self.rb_sinkers.toggled.connect(self.on_rope_mode_changed)
        self.rb_bottom_collar.toggled.connect(self.on_rope_mode_changed)

        # ✅ 초기 상태 반영
        self.on_rope_mode_changed()

        self.tab_collar.setLayout(form_collar)
        self.tab_widget.addTab(self.tab_collar, "Collar")



###########################################################################
#                              Mooring 탭                                 #
###########################################################################

        self.tab_mooring = QWidget()
        form_mooring = QFormLayout()

        self.mooring_sides_input = QLineEdit("4")
        self.mooring_circum_input = QLineEdit("200")
        self.mooring_side_nodes_input = QLineEdit("5")
        self.mooring_frame_z = QLineEdit("0")

        form_mooring.addRow(QLabel("Mooring Frame Sides"), self.mooring_sides_input)
        form_mooring.addRow(QLabel("Mooring Frame Circumference (m)"), self.mooring_circum_input)
        form_mooring.addRow(QLabel("Number of nodes per side"), self.mooring_side_nodes_input)
        form_mooring.addRow(QLabel("Z Mooring Frame (m)"), self.mooring_frame_z)

        moor_frame_btn = QPushButton("Submit Mooring Frame")
        moor_frame_btn.clicked.connect(self.submit_mooring_frame)
        form_mooring.addRow(moor_frame_btn)

        open_link_window_btn = QPushButton("Open Bridle line Generator")
        open_link_window_btn.clicked.connect(self.open_link_window)
        form_mooring.addRow(open_link_window_btn)

        hline1 = QFrame()
        hline1.setFrameShape(QFrame.HLine)
        hline1.setFrameShadow(QFrame.Sunken)
        form_mooring.addRow(hline1)

        open_mooring_rig_btn = QPushButton("Open Bouy & Anchor line")
        open_mooring_rig_btn.clicked.connect(self.open_mooring_rig_window)
        form_mooring.addRow(open_mooring_rig_btn)

        form_mooring.addRow(self.create_common_buttons())

        self.tab_mooring.setLayout(form_mooring)
        self.tab_widget.addTab(self.tab_mooring, "Mooring")

        frame.setLayout(h_layout)

    def open_side_rope_link_window(self):
        # bottom_collar(링) 데이터가 있어야 열 수 있음
        try:
            with open("Fish_Cage/bottom_collar_temp.json") as f:
                _ = json.load(f)
        except Exception as e:
            print("❌ bottom_collar_temp.json이 없습니다. 먼저 'Generate Bottom Collar'를 실행하세요.")
            return
        global BOOT_FLAG
        BOOT_FLAG = 1
        self.side_rope_link_window = UniversalLinkWindow(self)
        self.side_rope_link_window.show()

    def open_mooring_rig_window(self):
        try:
            # mooring이 있어야 진행 (없으면 경고)
            with open("Fish_Cage/mooring_temp.json", "r", encoding="utf-8") as f:
                _ = json.load(f)
        except Exception:
            QMessageBox.warning(self, "No mooring frame", "먼저 'Submit Mooring Frame'을 실행해 Mooring Frame을 만들어주세요.")
            return
        self.mooring_rig_window = MooringRigWindow(self)
        self.mooring_rig_window.show()

    def on_rope_mode_changed(self):
        """Sinker array면 circumference 비활성, Bottom collar면 활성."""
        use_bottom_collar = self.rb_bottom_collar.isChecked()
        self.bottomcollar_cir.setEnabled(use_bottom_collar)
        self.lbl_bottomcollar_cir.setEnabled(use_bottom_collar)
        self.bottomcollar_z.setEnabled(use_bottom_collar)
        self.lbl_bottomcollar_z.setEnabled(use_bottom_collar)
        self.btn_gen_bottom_collar.setEnabled(use_bottom_collar)
        self.open_side_rope_link_btn.setEnabled(use_bottom_collar)
        self.open_side_rope_link_sub_btn.setEnabled(use_bottom_collar)

    def open_link_window(self):
        # Bridle line Generator 버튼이 여기를 호출하므로 기본 프리셋 지정
        global BOOT_FLAG
        BOOT_FLAG = 3
        self.link_window = UniversalLinkWindow(self)
        self.link_window.show()

    def create_common_buttons(self):
        btn_widget = QWidget()
        vbox = QVBoxLayout()

        export_btn = QPushButton("Export Net to .med")
        export_btn.clicked.connect(self.export_and_exit)
        vbox.addWidget(export_btn)

        btn_widget.setLayout(vbox)
        return btn_widget

    def update_vtk_scene(self, points, edges):
        self.last_points = points  # ✅ 가장 최근 points 저장

        # 👇 안전 접근(존재할 때만 제거)
        actor = getattr(self, "actor", None)
        if actor:
            self.renderer.RemoveActor(actor)
        link_actor = getattr(self, "link_actor", None)
        if link_actor:
            self.renderer.RemoveActor(link_actor)

        vtk_points = vtk.vtkPoints()
        for p in points:
            vtk_points.InsertNextPoint(p)

        lines = vtk.vtkCellArray()
        for s, e in edges:
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, s)
            line.GetPointIds().SetId(1, e)
            lines.InsertNextCell(line)

        poly_data = vtk.vtkPolyData()
        poly_data.SetPoints(vtk_points)
        poly_data.SetLines(lines)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(poly_data)

        self.actor = vtk.vtkActor()
        self.actor.SetMapper(mapper)
        self.actor.GetProperty().SetColor(1, 0, 0)

        self.renderer.AddActor(self.actor)
        self.renderer.ResetCamera()
        self.vtkWidget.GetRenderWindow().Render()

    def submit_geometry(self):
        sides = int(self.sides_input.text())
        f_circum = float(self.f_circum_input.text())
        z_float = float(self.z_float_input.text())
        w_circum = float(self.bn_circum_input.text())
        z_weight = float(self.z_weight_input.text())

        points, edges = generate_cage_frame(sides, f_circum, z_float, w_circum, z_weight)
        self.update_vtk_scene(points, edges)

        data = {
            "sides": sides,
            "f_circumference": f_circum,
            "z_float": z_float,
            "w_circumference": w_circum,
            "z_weight": z_weight
        }
        with open("Fish_Cage/geometry_temp.json", "w") as f:
            json.dump(data, f)
        print("✅ geometry_temp.json saved!")

    def submit_net_parameter(self):
        with open("Fish_Cage/geometry_temp.json") as f:
            geo = json.load(f)

        half_mesh = float(self.half_mesh_input.text())
        diameter = float(self.diameter_input.text())
        bottomnetangle = float(self.bottomnetangle_input.text())

        points, edges, groups = generate_net_nodes_and_edges(
            geo["sides"],
            geo["f_circumference"],
            geo["w_circumference"],
            geo["z_float"],
            geo["z_weight"],
            half_mesh,
            bottomnetangle
        )
        self.update_vtk_scene(points, edges)


    def submit_mooring_frame(self):
        try:
            sides = int(self.mooring_sides_input.text())
            circumference = float(self.mooring_circum_input.text())
            nodes_per_side = int(self.mooring_side_nodes_input.text())
            z_level = float(self.mooring_frame_z.text())
            if nodes_per_side < 2:
                print("❌ Number of nodes per side must be ≥ 2")
                return
        except Exception:
            print("❌ Mooring Frame 입력 오류")
            return

        # ✅ 인자 순서: (sides, circumference, z_height, nodes_per_side)
        m_points, m_edges = generate_mooring_frame(
            sides, circumference, z_level, nodes_per_side
        )

        if not hasattr(self, "scene_dict"):
            self.scene_dict = {}
        self.scene_dict["mooring_frame"] = {
            "points": m_points,
            "edges": m_edges,
            "name": "MOORING_FRAME"
        }

        # ✅ 파일 저장 (변수명 통일)
        import os, json
        os.makedirs("Fish_Cage", exist_ok=True)
        with open("Fish_Cage/mooring_temp.json", "w", encoding="utf-8") as f:
            json.dump({"points": m_points, "edges": m_edges, "z_frame": z_level},
                    f, indent=2, ensure_ascii=False)
        print("✅ mooring_temp.json saved.")

        # ✅ 즉시 리빌드
        self.rebuild_main_scene()


    def submit_mooring_floater_link(self, points):
        if hasattr(self, 'link_actor') and self.link_actor:
            self.renderer.RemoveActor(self.link_actor)

        try:
            with open("Fish_Cage/link_temp.json") as f:
                link = json.load(f)
            link_edges = link.get("links", [])
        except:
            link_edges = []

        vtk_points = vtk.vtkPoints()
        for p in points:
            vtk_points.InsertNextPoint(p)

        lines = vtk.vtkCellArray()
        for s, e in link_edges:
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, s)
            line.GetPointIds().SetId(1, e)
            lines.InsertNextCell(line)

        poly = vtk.vtkPolyData()
        poly.SetPoints(vtk_points)
        poly.SetLines(lines)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(poly)

        self.link_actor = vtk.vtkActor()
        self.link_actor.SetMapper(mapper)
        self.link_actor.GetProperty().SetColor(1, 0, 0)

        self.renderer.AddActor(self.link_actor)
        self.vtkWidget.GetRenderWindow().Render()

    # gui_interface.py
    def export_and_exit(self):
        # 1) VTU는 GUI에서 즉시 저장
        try:
            out_dir = os.path.join("Fish_Cage", "saves")
            os.makedirs(out_dir, exist_ok=True)
            vtu_path = os.path.join(out_dir, "fish_cage.vtu")
            export_all_to_vtu(vtu_path)
            print(f"✅ VTU saved to: {os.path.abspath(vtu_path)}")
        except Exception as e:
            print(f"❌ VTU export failed: {e}")

        # 2) MED는 기존처럼 Salome Shell에서 export_med.py 실행
        print("👉 Now run export_med.py in Salome Shell to generate the MED file.")
        self.close()



    def open_side_rope_link_sub_window(self):
        # bottom_collar(링) 데이터가 있어야 열 수 있음
        try:
            with open("Fish_Cage/bottom_collar_temp.json") as f:
                _ = json.load(f)
        except Exception as e:
            print("❌ bottom_collar_temp.json이 없습니다. 먼저 'Generate Bottom Collar'를 실행하세요.")
            return
        # 부트 플래그 2: SIDE_ROPE_SUB 프리셋
        global BOOT_FLAG
        BOOT_FLAG = 2
        self.side_rope_link_sub_window = UniversalLinkWindow(self)
        self.side_rope_link_sub_window.show()

    def generate_outer_floater_json(self):
        try:
            sides            = int(self.outside_sides.text())
            bracket_per_side = int(self.bracket_num.text())
            circumference    = float(self.float_outside_cir.text())
            z_pos            = float(self.z_float_input.text())
        except ValueError:
            print("❌ 입력 오류")
            return

        # Net 로드
        try:
            with open("Fish_Cage/net_points.json") as f:
                net = json.load(f)
            net_pts, net_edges = net["points"], net["edges"]
        except:
            print("❌ net_points.json 불러오기 실패")
            return

        # 외부 floater 생성 (sides + bracket_per_side 기반)
        outer_pts, outer_edges = generate_outer_floater(
            sides=sides,
            bracket_per_side=bracket_per_side,
            circumference=circumference,
            z_pos=z_pos
        )

        # 병합 + VTK 업데이트
        offset     = len(net_pts)
        all_pts    = net_pts + outer_pts
        all_edges  = net_edges + [(s + offset, e + offset) for s, e in outer_edges]
        self.update_vtk_scene(all_pts, all_edges)

        print("✅ float_temp.json saved!")

    def open_floating_link_window(self):
        self.floating_link_window = UniversalLinkWindow(self)
        self.floating_link_window.show()   
        
    def submit_floating_floater_link(self):
        self.rebuild_main_scene()

    def generate_bottom_collar(self):
        """
        의도:
        - (우선) Net 에 '동일한 z' + '동일한 둘레(반지름)'를 가진 링이 충분히 있으면 그것을 bottom collar 로 인식(스냅)
        - (그렇지 않으면) 새 링을 생성하여 bottom collar 로 만든다
        - 한 번 클릭으로 저장+표시까지 완료
        """
        try:
            import math, json, os, numpy as np

            # ---- 입력 ----
            N      = int(self.outside_sides.text())              # 바닥 링의 분할(Outside Floater Sides)
            C_ring = float(self.bottomcollar_cir.text())
            z_bot  = float(self.bottomcollar_z.text())
            if N < 3:
                print("❌ Outside Floater Sides는 3 이상")
                return
            R = C_ring / (2 * np.pi)

            # ---- 데이터 로드 ----
            with open("Fish_Cage/net_points.json") as f:
                net = json.load(f)
            net_pts  = net.get("points", [])
            net_edges = net.get("edges", [])

            try:
                with open("Fish_Cage/float_temp.json") as f:
                    flo = json.load(f)
                float_pts = flo.get("points", [])
            except Exception:
                float_pts = []

            # ---- 오프셋 ----
            off_net   = len(net_pts)
            off_float = off_net + len(float_pts)

            # ---- 토러런스 ----
            EPS_Z      = 1e-6                    # z 일치
            EPS_R      = max(1e-3 * R, 1e-6)     # 반지름 일치(상대식)
            EPS_CENTER = max(1e-3 * R, 1e-5)     # 중심 배제

            def radius_xy(p): return math.hypot(p[0], p[1])

            # ---- 후보: z≈z_bot 인 '기존 Net' 노드(반지름 제한 X) ----
            z_plane_idx = []
            angles = []
            radii = []
            for i, p in enumerate(net_pts):
                if abs(p[2] - z_bot) > EPS_Z:
                    continue
                r = radius_xy(p)
                if r < EPS_CENTER:
                    continue
                z_plane_idx.append(i)
                angles.append(math.atan2(p[1], p[0]) % (2*math.pi))
                radii.append(r)

            # z 평면에 노드가 너무 적으면 신규 생성
            if len(z_plane_idx) < 3:
                # ▶ 신규 생성(N 분할)
                ring_global_idx = []
                points_new = []
                for k in range(N):
                    ang = 2*np.pi*k/N
                    x, y = R*np.cos(ang), R*np.sin(ang)
                    points_new.append([x, y, z_bot])
                    ring_global_idx.append(off_float + len(points_new) - 1)
                snap_mode = "NEW"
            else:
                # ▶ 각도 구간별 '최대 반지름' 노드 선택 → 외곽 링 복원
                #    (빈 수를 크게 잡아 변 사이 노드까지 포함되도록 한다)
                M_bins = max(180, min(8*len(z_plane_idx), 2000))
                bin_w = (2*math.pi) / M_bins

                # 각도 -> (최대r, idx) 매핑
                best_by_bin = [(-1.0, -1) for _ in range(M_bins)]
                for idx, (ang, r) in enumerate(zip(angles, radii)):
                    b = int(ang // bin_w)
                    if b >= M_bins:
                        b = M_bins - 1
                    if r > best_by_bin[b][0]:
                        best_by_bin[b] = (r, z_plane_idx[idx])

                # 유효한 후보만 추출
                ring_candidates = [idx for (rmax, idx) in best_by_bin if idx >= 0]
                # 중복 제거(안정적으로 각도 기준 정렬)
                cand_unique = sorted(set(ring_candidates), key=lambda i: math.atan2(net_pts[i][1], net_pts[i][0]) % (2*math.pi))

                # (선택) 목표 둘레와의 일치성 검사 → 원치 않으면 생략 가능
                #   R_eff = 중앙값 반지름, C_eff = 2πR_eff
                if C_ring > 0:
                    r_eff = np.median([radius_xy(net_pts[i]) for i in cand_unique])
                    c_eff = 2*math.pi*r_eff
                    if abs(c_eff - C_ring) > max(0.02*C_ring, 1e-6):
                        # 둘레가 너무 다르면 신규 생성으로 대체
                        ring_global_idx = []
                        points_new = []
                        for k in range(N):
                            ang = 2*np.pi*k/N
                            x, y = R*np.cos(ang), R*np.sin(ang)
                            points_new.append([x, y, z_bot])
                            ring_global_idx.append(off_float + len(points_new) - 1)
                        snap_mode = "NEW"
                    else:
                        ring_global_idx = cand_unique
                        points_new = []
                        snap_mode = "SNAP_FULLZ"
                else:
                    # 둘레 비교 안 할 경우 그대로 채택
                    ring_global_idx = cand_unique
                    points_new = []
                    snap_mode = "SNAP_FULLZ"

            # ---- 링 엣지(순환 연결) ----
            ring_edges = []
            M = len(ring_global_idx)
            for k in range(M):
                a = ring_global_idx[k]
                b = ring_global_idx[(k+1) % M]
                if a != b:
                    ring_edges.append((a, b))

            # (참고) 원래 Net 엣지에서 동일 엣지를 빼고 소유권을 옮기고 싶다면 여기서 정리 가능
            # ring_edge_set = set(tuple(sorted(e)) for e in ring_edges if e[0] < off_net and e[1] < off_net)
            # net_edges_filtered = [e for e in net_edges if tuple(sorted(e)) not in ring_edge_set]

            # ---- 저장 ----
            with open("Fish_Cage/bottom_collar_temp.json", "w") as f:
                json.dump({
                    "points_new": points_new,          # 스냅이면 일부/전부 빈 리스트일 수 있음
                    "ring_indices": ring_global_idx,   # 글로벌 인덱스(기존+신규 혼합)
                    "ring_edges": ring_edges
                }, f, indent=2)
            print(f"✅ bottom_collar_temp.json saved (mode={'SNAP' if snap_mode else 'NEW'}).")

            # ---- 즉시 표시 ----
            self.rebuild_main_scene()

        except Exception as e:
            print(f"❌ generate_bottom_collar 오류: {e}")


    def submit_side_rope_link(self):
            self.rebuild_main_scene()

    def rebuild_main_scene(self):

        # ---------- 1) Net ----------
        with open("Fish_Cage/net_points.json") as f:
            net = json.load(f)
        net_pts  = net.get("points", [])
        net_edges = [tuple(e) for e in net.get("edges", [])]

        points = list(net_pts)
        edges  = list(net_edges)

        # ---------- 2) Floating collar ----------
        float_pts, float_edges = [], []
        try:
            with open("Fish_Cage/float_temp.json") as f:
                flo = json.load(f)
            float_pts   = flo.get("points", [])
            float_edges = [tuple(e) for e in flo.get("edges", [])]
            off = len(points)
            points += float_pts
            edges  += [(s+off, e+off) for (s,e) in float_edges]
        except Exception:
            pass

        # ---------- 3) Bottom collar ----------
        # 포인트는 항상 추가(유니버설 링크가 참조할 수 있으므로),
        # 링 엣지는 라디오가 체크된 경우에만 표시
        bc_points_new, bc_ring_edges = [], []
        try:
            with open("Fish_Cage/bottom_collar_temp.json") as f:
                bc = json.load(f)
            bc_points_new = bc.get("points_new", [])
            bc_ring_edges = [tuple(e) for e in bc.get("ring_edges", [])]
        except Exception:
            pass

        points += bc_points_new
        if getattr(self, "rb_bottom_collar", None) and self.rb_bottom_collar.isChecked():
            edges += bc_ring_edges


        # ---------- 4) Mooring frame ----------
        moor_pts, moor_edges = [], []
        try:
            with open("Fish_Cage/mooring_temp.json") as f:
                moor = json.load(f)
            moor_pts   = moor.get("points", [])
            moor_edges = [tuple(e) for e in moor.get("edges", [])]
            off = len(points)
            points += moor_pts
            edges  += [(s+off, e+off) for (s,e) in moor_edges]
        except Exception:
            moor_pts = []
            pass

        # ---------- (옵션) 4) Floating bracket / side-rope 링크(기존 파일 사용 시) ----------
        try:
            with open("Fish_Cage/floating_link_temp.json") as f:
                fl = json.load(f)
            edges += [tuple(e) for e in fl.get("links", [])]   # 이미 글로벌 인덱스
        except Exception:
            pass

        try:
            with open("Fish_Cage/side_rope_link_temp.json") as f:
                sr = json.load(f)
            edges += [tuple(e) for e in sr.get("links", [])]   # 이미 글로벌 인덱스
        except Exception:
            pass

        # ---------- 5) Universal link ----------

        try:
            base_for_shift = len(net_pts) + len(float_pts) + len(bc_points_new) + len(moor_pts)
            curr_baseN = base_for_shift
            for path in sorted(glob.glob(os.path.join(LINKS_DIR, "*.json"))):
                with open(path, "r", encoding="utf-8") as f:
                    data = json.load(f)
                extra = data.get("extra_points", [])
                eds   = [tuple(int(x) for x in e) for e in data.get("edges", [])]
                saved_baseN = int(data.get("baseN", 0))

                # 저장 당시 baseN과 현재 baseN 차이를 보정해서 엣지 인덱스 시프트
                shift = curr_baseN - saved_baseN if saved_baseN else 0
                if shift:
                    eds = [
                        (i + (shift if i >= saved_baseN else 0),
                        j + (shift if j >= saved_baseN else 0))
                        for (i, j) in eds
                    ]

                points += extra
                edges  += eds
                curr_baseN += len(extra)
            print("✅ linked groups loaded from Fish_Cage/links")
        except Exception as e:
            print(f"ℹ️ link-group load skipped: {e}")

                # ---------- 6) Mooring frame ----------
        try:
            with open("Fish_Cage/mooring_temp.json") as f:
                moor = json.load(f)
            moor_pts   = moor.get("points", [])
            moor_edges = [tuple(e) for e in moor.get("edges", [])]
            off = len(points)
            points += moor_pts
            edges  += [(s+off, e+off) for (s,e) in moor_edges]
        except Exception:
            pass

        # ---------- 7) Draw ----------
        self.update_vtk_scene(points, edges)
        print("✅ Rebuilt main scene (with link_universal_temp).")

    def open_floating_bracket_link(self):
        global BOOT_FLAG
        BOOT_FLAG = 0
        self.universal_link_window = UniversalLinkWindow(self)
        self.universal_link_window.show()

    # =========================
    # ▶ 카메라 축뷰 헬퍼
    # =========================
    def _view_axis(self, axis: str):
        """
        axis in {'x','y','z'}.
        Shift 키와 함께 누르면 반대방향(-X/-Y/-Z)으로 봅니다.
        """
        try:
            cam = self.renderer.GetActiveCamera()
        except Exception:
            return
        if cam is None:
            return
        fp = cam.GetFocalPoint()
        dist = cam.GetDistance() or 100.0

        # Shift 누르면 반대 방향
        sign = -1.0 if (QApplication.keyboardModifiers() & Qt.ShiftModifier) else 1.0
        x, y, z = fp
        axis = (axis or '').lower()
        if axis == 'x':
            pos = (x + sign * dist, y, z); up = (0, 0, 1)
        elif axis == 'y':
            pos = (x, y + sign * dist, z); up = (0, 0, 1)
        elif axis == 'z':
            pos = (x, y, z + sign * dist); up = (0, 1, 0)  # 위에서 내려다보는 뷰
        else:
            return
        cam.SetPosition(*pos); cam.SetFocalPoint(*fp); cam.SetViewUp(*up)
        self.renderer.ResetCameraClippingRange(); self.vtkWidget.GetRenderWindow().Render()

    # =========================
    # ▶ 오른쪽 패널 토글 버튼 설치
    # =========================
    def _install_right_panel_toggle(self, splitter: QSplitter, handle_index: int = 1, default_right: int = 320):
        """
        splitter: [왼쪽 위젯 | 오른쪽 위젯] 순서의 QSplitter
        handle_index: 왼쪽/오른쪽 사이 손잡이 인덱스 (2개면 1)
        default_right: 복원값이 없을 때 오른쪽 패널 기본 폭(px)
        """
        handle = splitter.handle(handle_index)

        btn = QPushButton(">>", handle)
        btn.setFixedSize(26, 26)
        btn.setCursor(Qt.PointingHandCursor)
        btn.raise_()

        # 상태 저장
        splitter._right_collapsed = False
        splitter._stored_sizes = None

        def _recenter():
            # 손잡이 중앙 배치
            x = (handle.width() - btn.width()) // 2
            y = (handle.height() - btn.height()) // 2
            btn.move(max(0, x), max(0, y))

        def _toggle():
            sizes = splitter.sizes()
            left_idx = handle_index - 1
            right_idx = handle_index

            if not splitter._right_collapsed:
                # 접기: 현재 사이즈 저장 후 오른쪽을 0으로
                splitter._stored_sizes = sizes[:]
                total = sum(sizes) if sum(sizes) > 0 else (default_right * 2)
                new = sizes[:]
                new[right_idx] = 0
                new[left_idx] = max(1, total - new[right_idx])
                splitter.setSizes(new)
                splitter._right_collapsed = True
                btn.setText("<<")  # 다시 펼치기 표시
            else:
                # 펼치기: 저장 사이즈가 있으면 복원, 없으면 기본폭
                if splitter._stored_sizes and sum(splitter._stored_sizes) > 0:
                    splitter.setSizes(splitter._stored_sizes)
                else:
                    total = sum(splitter.sizes())
                    total = total if total > 0 else (default_right * 2)
                    left = max(1, total - default_right)
                    splitter.setSizes([left, default_right])
                splitter._right_collapsed = False
                btn.setText(">>")

        btn.clicked.connect(_toggle)

        # 이벤트: 리사이즈/표시/스플리터 이동 시 중앙 유지
        handle.installEventFilter(self)
        if not hasattr(self, "_toggle_handles"):
            self._toggle_handles = []
        self._toggle_handles.append((handle, _recenter))

        splitter.splitterMoved.connect(lambda *_: _recenter())
        QTimer.singleShot(0, _recenter)

    # 이벤트 필터: 스플리터 손잡이 리사이즈 때 버튼 중앙 정렬 유지
    def eventFilter(self, obj, ev):
        if hasattr(self, "_toggle_handles"):
            for h, recenter in self._toggle_handles:
                if obj is h and ev.type() in (QEvent.Resize, QEvent.Show):
                    recenter()
        return super().eventFilter(obj, ev)
    
    def _install_axes_gizmo(self, interactor, viewport=(0.82, 0.82, 1.0, 1.0)):
        """
        VTK 좌하단에 XYZ 축 표시(축 길이/색상 기본값 사용).
        viewport: (xmin, ymin, xmax, ymax) 0~1 정규 좌표
        """
        axes = vtkAxesActor()  # 기본: X=Red, Y=Green, Z=Blue
        # 필요 시 크기 조절: axes.SetTotalLength(1.0, 1.0, 1.0)

        w = vtkOrientationMarkerWidget()
        w.SetOrientationMarker(axes)
        w.SetInteractor(interactor)
        w.SetViewport(*viewport)   # 좌하단 18% 박스
        w.EnabledOn()
        w.InteractiveOff()         # 뷰 조작 방해하지 않도록

        # 참조 유지(가비지 컬렉션 방지)
        self._axes_actor = axes
        self._axes_widget = w

    # ---------- Data reset helpers ----------
    def _msg_yes(self, title, text):
        ret = QMessageBox.question(self, title, text,
                                   QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        return ret == QMessageBox.Yes

    def _delete_confirm(self, paths, question_text):
        if not self._msg_yes("Confirm delete", question_text):
            return
        deleted, missing, failed = [], [], []
        import os
        for p in paths:
            try:
                if os.path.exists(p):
                    os.remove(p)
                    deleted.append(p)
                else:
                    missing.append(p)
            except Exception as e:
                failed.append((p, str(e)))
        print(f"🗑️ deleted: {deleted}")
        if missing:
            print(f"ℹ️ missing (skipped): {missing}")
        if failed:
            print(f"⚠️ failed: {failed}")
        self._refresh_after_data_change()

    def _delete_links_dir(self, links_dir):
        import os, glob
        if not self._msg_yes("Confirm delete", "links/*.json을 모두 삭제할까요?"):
            return
        if not os.path.isdir(links_dir):
            print("ℹ️ links 폴더가 없어 생략")
            self._refresh_after_data_change()
            return
        failed = []
        for f in glob.glob(os.path.join(links_dir, "*.json")):
            try:
                os.remove(f)
            except Exception as e:
                failed.append((f, str(e)))
        if failed:
            print(f"⚠️ 일부 링크 파일 삭제 실패: {failed}")
        self._refresh_after_data_change()

    def _full_reset_fish_cage(self):
        """
        Fish_Cage 디렉토리 내의 모든 JSON(하위 links 포함)을 삭제합니다.
        디렉토리는 유지하고, links 폴더는 비운 뒤 다시 생성합니다.
        """
        import os, glob, shutil
        if not self._msg_yes("✅ 전체 초기화",
                             "Fish_Cage 폴더 내의 모든 .json 파일을 삭제합니다. 진행할까요?"):
            return
        base = "Fish_Cage"
        if not os.path.isdir(base):
            print("ℹ️ Fish_Cage 폴더가 없어 생략.")
            self._refresh_after_data_change()
            return

        # 1) 최상위 JSON 삭제
        top_failed = []
        for f in glob.glob(os.path.join(base, "*.json")):
            try:
                os.remove(f)
            except Exception as e:
                top_failed.append((f, str(e)))

        # 2) links/*.json 비우기 (폴더는 유지)
        links_dir = os.path.join(base, "links")
        if os.path.isdir(links_dir):
            for f in glob.glob(os.path.join(links_dir, "*.json")):
                try:
                    os.remove(f)
                except Exception as e:
                    top_failed.append((f, str(e)))
        else:
            os.makedirs(links_dir, exist_ok=True)

        if top_failed:
            print(f"⚠️ 일부 파일 삭제 실패: {top_failed}")
        else:
            print("🧹 전체 초기화 완료: Fish_Cage 내 JSON 삭제")

        self._refresh_after_data_change()

    def _clear_scene(self):
        """VTK 씬 비우기(리빌드 실패 시 대비)."""
        try:
            if self.actor:
                self.renderer.RemoveActor(self.actor)
                self.actor = None
            if hasattr(self, "link_actor") and self.link_actor:
                self.renderer.RemoveActor(self.link_actor)
                self.link_actor = None
            self.vtkWidget.GetRenderWindow().Render()
        except Exception:
            pass

    def _refresh_after_data_change(self):
        """
        파일 삭제/초기화 이후 장면 리프레시.
        net_points.json이 없어도 죽지 않게 예외 처리.
        """
        try:
            self.rebuild_main_scene()
        except Exception as e:
            print(f"ℹ️ rebuild_main_scene 실패 → 씬 초기화: {e}")
            self._clear_scene()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
