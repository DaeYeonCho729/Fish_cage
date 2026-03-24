# simulation_gui_interface.py
import sys
import vtk
import os, json
from vtkmodules.util.numpy_support import vtk_to_numpy
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QSplitter, QTabWidget, QLabel,
    QFileDialog, QGroupBox, QComboBox, QDoubleSpinBox, QPushButton, QFormLayout, QMessageBox, QVBoxLayout, QCheckBox, QLineEdit
)
from PyQt5.QtCore import Qt, QLocale
from vtkmodules.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera
from PyQt5.QtGui import QDoubleValidator

EDGE_GROUP_ORDER = [
    "floater", "net twine", "bottom ring",
    "sinker",
    "side ropes", "bridle line",
    "buoyline1", "buoyline2",
    "anchor line 1", "anchor line 2",
    "mooring frame", "bracket",
]
EDGE_GROUP_LABEL = {
    "floater":"Floater", "net twine":"Net twine", "bottom ring":"Bottom ring",
    "sinker":"Sinker",
    "side ropes":"Side ropes", "bridle line":"Bridle line",
    "buoyline1":"Buoy line 1", "buoyline2":"Buoy line 2",
    "anchor line 1":"Anchor line A", "anchor line 2":"Anchor line B",
    "mooring frame":"Mooring frame", "bracket":"Bracket",
}

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.group_checks = {}

        self.setWindowTitle("Fish Cage Simulation GUI (DaeYeon Cho)")
        self.setMinimumSize(1200, 800)
        self._build_ui()

    def _build_ui(self):
        self._form_state = {}      # { "floater": {...}, "net twine": {...}, "bottom ring": {...} }
        self._current_group = None # 마지막 선택 그룹 키

        # ---- 중앙 위젯 + 레이아웃
        central = QWidget(self)
        layout = QHBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setCentralWidget(central)

        # ---- 스플리터
        self.splitter = QSplitter(Qt.Horizontal, central)
        layout.addWidget(self.splitter)

        # ---- (오른쪽 먼저) 탭 위젯 생성  ← ★ self.tab_widget를 먼저 만든다
        self.tab_widget = QTabWidget(self)

        # ↓↓↓ 여기 3줄 추가 (같은 들여쓰기: 8칸)
        run_tab = QWidget(self)
        self.tab_widget.addTab(run_tab, "Setting")
        self.attach_edge_material_form(QVBoxLayout(run_tab))

        #environment 탭 추가
        env_tab = QWidget(self)
        self.tab_widget.addTab(env_tab, "Environment")
        self.attach_environment_form(QVBoxLayout(env_tab))

        # ---- (왼쪽) VTK 뷰어 생성
        self.vtk_widget = QVTKRenderWindowInteractor(self)

        # --- ADD: MED/그룹/재료 상태 ---
        self.med_actor = None               # 현재 불러온 MED 메쉬 actor
        self.med_ugrid = None               # vtkUnstructuredGrid 보관
        self.edge_groups = {}               # { group_name: [cellIds] } (1D만)
        self.material_by_group = {}         # { group_name: {E:.., rho:.., d:..} }

        # --- /ADD ---

        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(0.1, 0.1, 0.12)

        # RenderWindow & Interactor
        self.vtk_widget.GetRenderWindow().AddRenderer(self.ren)
        interactor = self.vtk_widget.GetRenderWindow().GetInteractor()
        interactor.SetInteractorStyle(vtkInteractorStyleTrackballCamera())

        # ---- 스플리터에 위젯 추가 (좌: VTK, 우: 탭)
        self.splitter.addWidget(self.vtk_widget)
        self.splitter.addWidget(self.tab_widget)
        self.splitter.setSizes([900, 300])

        # ---- 테스트용 축 추가 (보이는지 확인)
        axes = vtk.vtkAxesActor()
        widget = vtk.vtkOrientationMarkerWidget()
        widget.SetOrientationMarker(axes)
        widget.SetViewport(0.0, 0.0, 0.15, 0.15)
        widget.SetInteractor(interactor)
        widget.EnabledOn()
        widget.InteractiveOff()

        self.ren.ResetCamera()
        self.vtk_widget.GetRenderWindow().Render()
        interactor.Initialize()

    def load_med_via_dialog(self):
        path, _ = QFileDialog.getOpenFileName(self, "Open MED file", "", "MED Files (*.med *.rmed)")
        if not path:
            return
        self.load_med_file(path)

    def load_vtu_via_dialog(self):
        from PyQt5.QtWidgets import QFileDialog, QMessageBox
        import vtk

        path, _ = QFileDialog.getOpenFileName(self, "Open VTU (PolyData)", "", "VTU Files (*.vtu)")
        if not path:
            return
        
        self.current_vtu_path = path

        # ---- PolyData 전용 리더 ----
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(path)
        try:
            reader.Update()
        except Exception as e:
            QMessageBox.warning(self, "VTU", f"읽기 실패: {e}")
            return

        poly = reader.GetOutput()
        if poly is None or poly.GetNumberOfCells() == 0:
            QMessageBox.warning(self, "VTU", "유효한 PolyData(.vtu)가 아닙니다.")
            return

        # 기존 actor 제거 (있을 때만)
        if hasattr(self, "med_actor") and self.med_actor is not None:
            try:
                self.ren.RemoveActor(self.med_actor)
            except Exception:
                pass
            self.med_actor = None

        # ---- 표시 ----
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(poly)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        self.ren.AddActor(actor)
        self.ren.ResetCamera()
        self.vtk_widget.GetRenderWindow().Render()

        # 상태 저장
        self.med_actor = actor
        self.med_poly  = poly   # ← PolyData로 저장

        # 그룹 업데이트 (edge_group / node_group)
        self.edge_groups, self.node_groups = self._extract_groups_from_polydata(poly)
        
        # edge_groups 만든 직후에 추가
        self.meshinfo = self._try_load_meshinfo(self.current_vtu_path)
        S_a = bool((self.meshinfo or {}).get("S_a", False))

        # S_a=True면 sinker 그룹을 bottom ring에서 파생
        if S_a and "bottom ring" in (self.edge_groups or {}) and self.edge_groups["bottom ring"]:
            self.edge_groups["sinker"] = list(self.edge_groups["bottom ring"])
        else:
            # 없으면 아예 비워두거나 키를 삭제 (UI enable을 위해 여기선 비움)
            self.edge_groups.pop("sinker", None)

        self._sync_group_ui()
        # UI 콤보가 있다면 갱신
        if hasattr(self, "cb_edge_group") and self.cb_edge_group is not None:
            self.cb_edge_group.clear()
            self.cb_edge_group.addItems(sorted(self.edge_groups.keys()))
        if hasattr(self, "cb_node_group") and self.cb_node_group is not None:
            self.cb_node_group.clear()
            self.cb_node_group.addItems(sorted(self.node_groups.keys()))

        self._set_net_section_enabled(False)
        self._set_floater_section_enabled(False)

    def _read_common_numbers(self):
        # 콤마 제거 후 float 캐스팅
        def _val(le): 
            t = (le.text() or "").replace(",", "").strip()
            return float(t) if t != "" else None
        return _val(self.in_E), _val(self.in_rho), _val(self.in_d)

    def _save_form_to_state(self, key: str):
        if not key: 
            return
        E, rho, d = self._read_common_numbers()
        if key == "floater":
            # floater는 d 미사용
            def _v(le):
                t = (le.text() or "").replace(",", "").strip()
                return float(t) if t != "" else None
            flo_diam = _v(self.in_floater_diam)
            flo_thk  = _v(self.in_floater_thk)
            self._form_state[key] = {
                "FloaterYoungModule": E, "RHO": rho,
                "Floaterdiameter": flo_diam, "Floaterthickness": flo_thk
            }
        elif key == "net twine":
            def _v(le):
                t = (le.text() or "").replace(",", "").strip()
                return float(t) if t != "" else None
            mesh_len = _v(self.in_mesh_length)
            bnw      = _v(self.in_bottom_net_weight)
            self._form_state[key] = {
                "NetYoungModule": E, "RHO": rho,
                "twine_dia": d, "Mesh_length": mesh_len, "BottomNet_cen_weight": bnw
            }
        elif key == "bottom ring":
            def _v(le):
                t = (le.text() or "").replace(",", "").strip()
                return float(t) if t != "" else None
            self._form_state[key] = {
                "BrYoungModule": E, "RHO": rho,
                "Br_diameter": d
            }
        elif key == "bracket":
            self._form_state[key] = {
                "BrktYoungModule": E,
                "RHO": rho,
                "Brkt_diameter": d,
            }
        elif key == "side ropes":
            self._form_state[key] = {
                "SideRopeYoungModule": E,
                "RHO": rho,
                "SideRope_diameter": d
            }
        elif key == "bridle line":
            self._form_state[key] = {
                "BridleYoungModule": E,
                "RHO": rho,
                "Bridle_diameter": d
            }
        elif key in ("buoyline1", "buoyline2"):
            def _v(le):
                t = (le.text() or "").replace(",", "").strip()
                return float(t) if t != "" else None

            bf = _v(getattr(self, "in_buoy_force", None))

            if key == "buoyline1":
                self._form_state[key] = {
                    "Buoy1YoungModule": E, "RHO": rho, "Buoy1_diameter": d,
                    "BuoyForce": bf
                }
            else:  # buoyline2
                self._form_state[key] = {
                    "Buoy2YoungModule": E, "RHO": rho, "Buoy2_diameter": d,
                    "BuoyForce": bf
                }
        elif key == "anchor line 1":
            self._form_state[key] = {
                "AnchorAYoungModule": E,
                "RHO": rho,
            "AnchorA_diameter": d
            }
        elif key == "anchor line 2":
            self._form_state[key] = {
                "AnchorBYoungModule": E,
                "RHO": rho,
                "AnchorB_diameter": d
            }
        elif key == "mooring line":
            self._form_state[key] = {
                "MooringYoungModule": E,
                "RHO": rho,
                "Mooring_diameter": d
            }
        elif key == "mooring frame":
            self._form_state[key] = {
                "MoorFrameYoungModule": E,
                "RHO": rho,
                "MoorFrame_diameter": d
            }


    def _load_state_to_form(self, key: str):
        # 먼저 전체 폼을 비우지 않습니다(입력 유지) — 상태에 값이 있으면 그 값으로 덮어쓰기
        st = self._form_state.get(key, {})
        def _set(le, v):
            if v is None: return
            # 콤마 포맷 재사용
            s = f"{v}"
            le.blockSignals(True); le.setText(s); le.blockSignals(False)
            self._format_with_comma(le, le.text())
        if key == "bracket":  # ✅ 추가 분기(키 이름이 다르니 직접 꺼내기)
            _set(self.in_E,   st.get("BrktYoungModule"))
            _set(self.in_rho, st.get("RHO"))
            _set(self.in_d,   st.get("Brkt_diameter"))
            return
        # 공통
        if "FloaterYoungModule" in st or "NetYoungModule" in st or "BrYoungModule" in st or "RHO" in st:
            # 어떤 섹션이든 E/rho 키가 있으면 공통칸 세팅
            _E = st.get("FloaterYoungModule", st.get("NetYoungModule", st.get("BrYoungModule")))
            _set(self.in_E, _E)
            _set(self.in_rho, st.get("RHO"))
        # 섹션별
        if key == "floater":
            _set(self.in_floater_diam, st.get("Floaterdiameter"))
            _set(self.in_floater_thk,  st.get("Floaterthickness"))
        elif key == "net twine":
            _set(self.in_d,                st.get("twine_dia"))
            _set(self.in_mesh_length,      st.get("Mesh_length"))
            _set(self.in_bottom_net_weight,st.get("BottomNet_cen_weight"))
        elif key == "bottom ring":
            _set(self.in_d,              st.get("Br_diameter"))

    def _extract_groups_from_polydata(self, poly):
        """
        export_vtu.py가 기록한:
        - CellData:  edge_group (vtkStringArray)
        - PointData: node_group (vtkStringArray)
        를 읽어 dict로 반환.
        """
        import vtk
        edge_groups = {}
        node_groups = {}

        # --- CellData → edge_group ---
        cdata = poly.GetCellData()
        if cdata is not None:
            arr = cdata.GetAbstractArray("edge_group")
            if arr is not None and isinstance(arr, vtk.vtkStringArray):
                ncell = poly.GetNumberOfCells()
                for ci in range(ncell):
                    name = arr.GetValue(ci) if ci < arr.GetNumberOfValues() else ""
                    name = name or ""
                    edge_groups.setdefault(name, []).append(ci)

        # --- PointData → node_group ---
        pdata = poly.GetPointData()
        if pdata is not None:
            arr = pdata.GetAbstractArray("node_group")
            if arr is not None and isinstance(arr, vtk.vtkStringArray):
                npt = poly.GetNumberOfPoints()
                for pi in range(npt):
                    name = arr.GetValue(pi) if pi < arr.GetNumberOfValues() else ""
                    name = name or ""
                    node_groups.setdefault(name, []).append(pi)

        return edge_groups, node_groups


    def get_loaded_polydata(self):
        """필요 시 다른 함수에서 현재 PolyData를 얻을 때 사용."""
        return getattr(self, "med_poly", None)


    def _debug_list_arrays(self, ugrid):
        cdat = ugrid.GetCellData()
        fdat = ugrid.GetFieldData()
        print("== CellData arrays ==")
        for i in range(cdat.GetNumberOfArrays()):
            arr = cdat.GetAbstractArray(i)
            print(f"  [{i}] name='{arr.GetName()}', type={type(arr).__name__}, N={arr.GetNumberOfTuples()}")
            if hasattr(arr, "GetValue"):
                # 문자열 배열이면 앞 10개 샘플
                try:
                    vals = [arr.GetValue(k) for k in range(min(10, arr.GetNumberOfTuples()))]
                    print("       sample:", vals)
                except Exception:
                    pass
        print("== FieldData arrays ==")
        for i in range(fdat.GetNumberOfArrays()):
            arr = fdat.GetAbstractArray(i)
            print(f"  [{i}] name='{arr.GetName()}', type={type(arr).__name__}, N={arr.GetNumberOfTuples()}")
            if hasattr(arr, "GetValue"):
                try:
                    vals = [arr.GetValue(k) for k in range(min(10, arr.GetNumberOfTuples()))]
                    print("       sample:", vals)
                except Exception:
                    pass
                
    def _build_vtk_from_meshio(self, path: str) -> vtk.vtkUnstructuredGrid:
        import meshio, numpy as np

        m = meshio.read(path)

        # 1) points
        pts = vtk.vtkPoints()
        pts.SetNumberOfPoints(len(m.points))
        for i, p in enumerate(m.points):
            x = float(p[0]); y = float(p[1]); z = float(p[2]) if len(p) > 2 else 0.0
            pts.SetPoint(i, x, y, z)

        ugrid = vtk.vtkUnstructuredGrid()
        ugrid.SetPoints(pts)

        # 2) cells + 블록 로컬→VTK 글로벌 인덱스 매핑
        cell_type_map = {
            "line":  vtk.VTK_LINE,
            "line3": vtk.VTK_QUADRATIC_EDGE,
            "triangle": vtk.VTK_TRIANGLE,
            "quad":     vtk.VTK_QUAD,
            "tetra":    vtk.VTK_TETRA,
            "hexahedron": vtk.VTK_HEXAHEDRON,
        }
        block_local_to_vtk = {}  # {"line": [vtkCellId,...], ...}
        for block in m.cells:
            ctype = block.type
            vtk_type = cell_type_map.get(ctype, None)
            if vtk_type is None:
                continue
            for conn in block.data:
                ids = vtk.vtkIdList()
                for nid in conn:
                    ids.InsertNextId(int(nid))
                vtk_id = ugrid.InsertNextCell(vtk_type, ids)
                block_local_to_vtk.setdefault(ctype, []).append(vtk_id)

        # 3) Group 문자열 배열 생성 (초기화: 빈 문자열)
        group_arr = vtk.vtkStringArray()
        group_arr.SetName("Group")
        group_arr.SetNumberOfValues(ugrid.GetNumberOfCells())
        for ci in range(ugrid.GetNumberOfCells()):
            group_arr.SetValue(ci, "")

        # --------- (A) cell_sets에서 그룹명 적용 ---------
        cell_sets = getattr(m, "cell_sets_dict", None) or getattr(m, "cell_sets", {})
        if isinstance(cell_sets, dict) and cell_sets:
            for gname, mapping in cell_sets.items():
                if isinstance(mapping, dict):  # 신형: {"line":[idx...], ...}
                    for ctype, idxs in mapping.items():
                        vtk_ids = block_local_to_vtk.get(ctype, [])
                        for local_i in idxs:
                            li = int(local_i)
                            if 0 <= li < len(vtk_ids):
                                group_arr.SetValue(int(vtk_ids[li]), str(gname))
                else:
                    # 구형: list/ndarray → 라인 타입에 일괄 적용
                    for ctype in ("line", "line3"):
                        vtk_ids = block_local_to_vtk.get(ctype, [])
                        for local_i in mapping:
                            li = int(local_i)
                            if 0 <= li < len(vtk_ids):
                                group_arr.SetValue(int(vtk_ids[li]), str(gname))

        # --------- (B) cell_data_dict에서 보강 (group/family류 키) ---------
        cd = getattr(m, "cell_data_dict", None)
        if isinstance(cd, dict):
            for key, mapping in cd.items():
                if str(key).lower() in ("group", "groups", "family", "familyname", "name"):
                    for ctype, values in mapping.items():
                        vtk_ids = block_local_to_vtk.get(ctype, [])
                        for local_i, val in enumerate(values):
                            if 0 <= local_i < len(vtk_ids):
                                group_arr.SetValue(int(vtk_ids[local_i]), str(val))

        # --------- (C) field_data에 그룹 메타가 있으면 매핑(일부 MED→meshio 케이스) ---------
        # meshio의 field_data는 {"GroupName": np.array([...meta...])} 형태일 수 있음
        # cell_sets가 비어 있고 field_data만 있는 경우를 대비하여,
        # 이름만이라도 유지되도록 (기존 값이 빈 문자열일 때) 남은 라인 셀에 균등 분배(최후수단).
        if hasattr(m, "field_data") and isinstance(m.field_data, dict):
            # 아직 이름이 안 들어간 라인 셀 추출
            unset_line_cells = []
            for ctype in ("line", "line3"):
                for vtk_id in block_local_to_vtk.get(ctype, []):
                    if group_arr.GetValue(int(vtk_id)) == "":
                        unset_line_cells.append(int(vtk_id))
            if unset_line_cells:
                # field_data 키들을 순회하며 이름만이라도 채운다 (균등 분배)
                names = list(m.field_data.keys())
                if names:
                    for idx, vtk_id in enumerate(unset_line_cells):
                        group_arr.SetValue(vtk_id, str(names[idx % len(names)]))

        ugrid.GetCellData().AddArray(group_arr)
        return ugrid


    def _try_load_meshinfo(self, vtu_path: str) -> dict:
        """
        vtu 경로 기준으로 meshinfo.py를 찾아서 meshinfo dict를 로드.
        (찾는 순서)
        1) vtu 폴더: meshinfo.py
        2) 상위 폴더: asterinput/meshinfo.py
        3) 더 상위로 계속 탐색
        """
        def _candidates(cur: str):
            return [
                os.path.join(cur, "meshinfo.py"),
                os.path.join(cur, "asterinput", "meshinfo.py"),
            ]

        cur = os.path.abspath(os.path.dirname(vtu_path))
        for _ in range(10):
            for p in _candidates(cur):
                if os.path.isfile(p):
                    scope = {}
                    try:
                        with open(p, "r", encoding="utf-8") as f:
                            code = f.read()
                        exec(code, scope)
                        mi = scope.get("meshinfo", {})
                        return mi if isinstance(mi, dict) else {}
                    except Exception:
                        return {}
            parent = os.path.dirname(cur)
            if parent == cur:
                break
            cur = parent
        return {}

    def _find_edge_groups(self, ugrid: vtk.vtkUnstructuredGrid):
        fixed_names = {
            "floater", "net twine", "bottom ring",
            "side ropes", "bridle line",
            "buoyline1", "buoyline2",
            "anchor line 1", "anchor line 2",
            "mooring frame", "bracket"
        }

        # 1D 셀만 선별
        edge_ids = []
        for cid in range(ugrid.GetNumberOfCells()):
            ctype = ugrid.GetCellType(cid)
            if ctype in (vtk.VTK_LINE, vtk.VTK_QUADRATIC_EDGE, vtk.VTK_CUBIC_LINE):
                edge_ids.append(cid)
        if not edge_ids:
            return {}

        cdat = ugrid.GetCellData()

        # CellData의 문자열 배열 "Group"/"Groups" 우선
        group_arr = None
        for i in range(cdat.GetNumberOfArrays()):
            arr = cdat.GetAbstractArray(i)
            if isinstance(arr, vtk.vtkStringArray):
                nm = (arr.GetName() or "").lower()
                if nm in ("group", "groups"):
                    group_arr = arr
                    break

        if group_arr is None:
            # 그룹 배열이 없을 때는 전체 1D를 묶어서라도 반환
            return {"EDGES_ALL": edge_ids}

        group_map = {}
        for cid in edge_ids:
            g = group_arr.GetValue(cid)
            if g in fixed_names:
                group_map.setdefault(g, []).append(cid)

        return group_map if group_map else {"EDGES_ALL": edge_ids}

    def _count_unique_points_from_cell_ids(self, data_obj, cell_ids):
        pts = set()
        idlist = vtk.vtkIdList()
        for cid in cell_ids:
            data_obj.GetCell(int(cid)).GetPointIds(idlist)
            for i in range(idlist.GetNumberOfIds()):
                pts.add(idlist.GetId(i))
        return len(pts)

    def _sync_group_ui(self):
        have = set((self.edge_groups or {}).keys())
        for key, cb in getattr(self, "group_checks", {}).items():
            if key in have and len(self.edge_groups.get(key, [])) > 0:  
                n = len(self.edge_groups[key])

                if key == "sinker" and getattr(self, "med_data", None) is not None:
                    try:
                        n = self._count_unique_points_from_cell_ids(self.med_data, self.edge_groups[key])
                    except Exception:
                        n = len(self.edge_groups[key])

                base = EDGE_GROUP_LABEL.get(key, key)
                cb.blockSignals(True)
                cb.setEnabled(True)
                cb.setText(f"{base}  ({n})")
                cb.blockSignals(False)
            else:
                base = EDGE_GROUP_LABEL.get(key, key)
                cb.blockSignals(True)
                cb.setChecked(False)
                cb.setEnabled(False)
                cb.setText(base)
                cb.blockSignals(False)

        # ===== sinker가 있으면 bottom ring 비활성화 =====
        sinker_cb = getattr(self, "group_checks", {}).get("sinker", None)
        bottom_cb = getattr(self, "group_checks", {}).get("bottom ring", None)
        if sinker_cb and bottom_cb and sinker_cb.isEnabled():
            bottom_cb.blockSignals(True)
            bottom_cb.setChecked(False)
            bottom_cb.setEnabled(False)
            bottom_cb.setText(EDGE_GROUP_LABEL.get("bottom ring", "bottom ring"))
            bottom_cb.blockSignals(False)

    def _on_group_toggled(self, key: str, checked: bool):
        if not hasattr(self, "group_checks"):
            return

        # === (A) 전환 직전, 이전 그룹 상태 저장 ===
        if checked:
            if self._current_group and self._current_group != key:
                self._save_form_to_state(self._current_group)

        # 1) 단일 선택 강제
        if checked:
            for k, cb in self.group_checks.items():
                if k != key and cb.isChecked():
                    cb.blockSignals(True)
                    cb.setChecked(False)
                    cb.blockSignals(False)

        # 2) 시각화
        if checked:
            self.preview_edge_group(key)
        else:
            if all(not cb.isChecked() for cb in self.group_checks.values()):
                self._clear_edge_preview()

        # 3) 섹션 enable
        net_on     = (key == "net twine" and checked)
        floater_on = (key == "floater" and checked)
        sinker_on  = (key == "sinker" and checked) 

        self._set_net_section_enabled(net_on)
        self._set_floater_section_enabled(floater_on)

        # grp_bottom(SINKER 박스)은 sinker일 때만 켜지게
        if hasattr(self, "grp_bottom"):
            self.grp_bottom.setEnabled(bool(sinker_on))

        # 입력칸도 sinker일 때만 켜지게
        w = getattr(self, "in_single_sinker_kg", None)
        if w is not None:
            w.setEnabled(bool(sinker_on))

        # 3-1) Floater 선택 시 d 비활성
        if hasattr(self, "in_d"):
            self.in_d.setEnabled(not floater_on)
            self.in_d.setPlaceholderText("" if floater_on else "0")

        # === (B) 새 그룹 로드 ===
        if checked:
            self._current_group = key
            self._load_state_to_form(key)
        else:
            if self._current_group == key:
                self._current_group = None

        if hasattr(self, "grp_buoy") and hasattr(self, "in_buoy_force"):
            on = bool(checked and key in ("buoyline1", "buoyline2"))
            self.grp_buoy.setEnabled(on)
            self.in_buoy_force.setEnabled(on)  # 선택: 명시적으로 같이
            self._buoy_force_target = key if on else None

    def _clear_edge_preview(self):
        if hasattr(self, "_edge_preview_actor") and self._edge_preview_actor:
            try:
                self.ren.RemoveActor(self._edge_preview_actor)
            except Exception:
                pass
            self._edge_preview_actor = None
            
        if hasattr(self, "_node_preview_actor") and self._node_preview_actor:
            try:
                self.ren.RemoveActor(self._node_preview_actor)
            except Exception:
                pass
            self._node_preview_actor = None

        if getattr(self, "med_actor", None):
            self.med_actor.GetProperty().SetColor(1,1,1)
            self.med_actor.GetProperty().SetOpacity(1.0)
        self.ren.Render()
        self.vtk_widget.GetRenderWindow().Render()

    def _set_bottom_section_enabled(self, on: bool):
        if hasattr(self, "grp_bottom"):
            self.grp_bottom.setEnabled(bool(False))
    
    def preview_edge_group(self, group_name: str):

        # ===== 모든 프리뷰(선/점) 먼저 제거 =====
        if hasattr(self, "_edge_preview_actor") and self._edge_preview_actor:
            try:
                self.ren.RemoveActor(self._edge_preview_actor)
            except Exception:
                pass
            self._edge_preview_actor = None

        if hasattr(self, "_node_preview_actor") and self._node_preview_actor:
            try:
                self.ren.RemoveActor(self._node_preview_actor)
            except Exception:
                pass
            self._node_preview_actor = None

        # 유효성 검사
        if not group_name or group_name not in getattr(self, "edge_groups", {}):
            return

        # 데이터 소스 (우선순위: PolyData → UnstructuredGrid)
        data_obj = getattr(self, "med_poly", None)
        if data_obj is None:
            data_obj = getattr(self, "med_ugrid", None)
        if data_obj is None:
            return

        ids = self.edge_groups[group_name]
        if not ids:
            return

        #sinker는 선이 아니라 "노드(점)"로 표시
        if group_name == "sinker":
            # 기존 하이라이트(선/점) 제거
            if hasattr(self, "_edge_preview_actor") and self._edge_preview_actor:
                try:
                    self.ren.RemoveActor(self._edge_preview_actor)
                except Exception:
                    pass
                self._edge_preview_actor = None
            if hasattr(self, "_node_preview_actor") and self._node_preview_actor:
                try:
                    self.ren.RemoveActor(self._node_preview_actor)
                except Exception:
                    pass
                self._node_preview_actor = None

            # 선택한 라인 셀들을 extract
            idarr = vtk.vtkIdTypeArray()
            idarr.SetNumberOfComponents(1)
            for cid in ids:
                idarr.InsertNextValue(int(cid))

            sel_node = vtk.vtkSelectionNode()
            sel_node.SetFieldType(vtk.vtkSelectionNode.CELL)
            sel_node.SetContentType(vtk.vtkSelectionNode.INDICES)
            sel_node.SetSelectionList(idarr)

            selection = vtk.vtkSelection()
            selection.AddNode(sel_node)

            extract = vtk.vtkExtractSelection()
            extract.SetInputData(0, data_obj)
            extract.SetInputData(1, selection)
            extract.Update()

            # extract 결과에서 "포인트만" glyph로 찍기
            geom = vtk.vtkGeometryFilter()
            geom.SetInputConnection(extract.GetOutputPort())
            geom.Update()

            poly = geom.GetOutput()
            vg = vtk.vtkVertexGlyphFilter()
            vg.SetInputData(poly)
            vg.Update()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(vg.GetOutputPort())

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(1.0, 0.0, 0.0)  # 🔴 빨강
            actor.GetProperty().SetPointSize(10)         # 필요하면 조절

            # 배경 메쉬 디밍(원하면 유지)
            if getattr(self, "med_actor", None):
                self.med_actor.GetProperty().SetColor(0.75, 0.75, 0.75)
                self.med_actor.GetProperty().SetOpacity(0.35)

            self._node_preview_actor = actor
            self.ren.AddActor(actor)
            self.ren.Render()
            self.vtk_widget.GetRenderWindow().Render()
            return

        # 이전 하이라이트 제거
        if hasattr(self, "_edge_preview_actor") and self._edge_preview_actor:
            try:
                self.ren.RemoveActor(self._edge_preview_actor)
            except Exception:
                pass
            self._edge_preview_actor = None

        # 선택 ID 배열
        idarr = vtk.vtkIdTypeArray()
        idarr.SetNumberOfComponents(1)
        for cid in ids:
            idarr.InsertNextValue(int(cid))

        # Selection → Extract
        sel_node = vtk.vtkSelectionNode()
        sel_node.SetFieldType(vtk.vtkSelectionNode.CELL)           # 셀 선택
        sel_node.SetContentType(vtk.vtkSelectionNode.INDICES)      # 인덱스 기반
        sel_node.SetSelectionList(idarr)

        selection = vtk.vtkSelection()
        selection.AddNode(sel_node)

        extract = vtk.vtkExtractSelection()
        extract.SetInputData(0, data_obj)
        extract.SetInputData(1, selection)
        extract.Update()

        # 표면/라인으로 변환해서 그리기 (데이터 타입 상관없이)
        geom = vtk.vtkGeometryFilter()
        geom.SetInputConnection(extract.GetOutputPort())
        geom.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(geom.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1.0, 0.0, 0.0)   # 🔴 빨간색
        actor.GetProperty().SetLineWidth(4.0)
        actor.GetProperty().SetOpacity(1.0)

        # 배경 메쉬 살짝 디밍 (원치 않으면 아래 2줄 주석 처리)
        if getattr(self, "med_actor", None):
            self.med_actor.GetProperty().SetColor(0.75, 0.75, 0.75)
            self.med_actor.GetProperty().SetOpacity(0.35)

        self._edge_preview_actor = actor
        self.ren.AddActor(actor)
        self.ren.Render()
        self.vtk_widget.GetRenderWindow().Render()


    def assign_material_to_group(self, group_name: str, E: float, rho: float, d: float):
        if not group_name:
            return
        self.material_by_group[group_name] = {"E": E, "rho": rho, "d": d}
        out = os.path.join(os.path.dirname(__file__), "edge_materials.json")
        try:
            with open(out, "w", encoding="utf-8") as f:
                json.dump(self.material_by_group, f, ensure_ascii=False, indent=2)
        except Exception as e:
            print("material save failed:", e)

    def attach_edge_material_form(self, parent_layout):
        w = QWidget(self)
        form = QFormLayout(w); form.setContentsMargins(6,6,6,6)

        btn_load = QPushButton("Load .vtu", w)
        btn_load.clicked.connect(self.load_vtu_via_dialog)
        form.addRow("VTU", btn_load)

        grp_box = QGroupBox("Edge groups")
        grp_lay = QVBoxLayout(grp_box)
        self.group_checks = {}

        for key in EDGE_GROUP_ORDER:
            cb = QCheckBox(EDGE_GROUP_LABEL.get(key, key))
            cb.setEnabled(False)
            cb.setChecked(False)
            cb.toggled.connect(lambda on, k=key: self._on_group_toggled(k, on))
            grp_lay.addWidget(cb)
            self.group_checks[key] = cb

        form.addRow(grp_box)

        # self.grp_env = QGroupBox("ENVIRONMENT")
        # self.grp_env.setCheckable(True)
        # self.grp_env.setChecked(False)
        # env_form = QFormLayout(self.grp_env)

        # self.in_dt = QLineEdit();        self.in_dt.setPlaceholderText("(s)")
        # self.in_duration = QLineEdit();  self.in_duration.setPlaceholderText("(s)")
        # self.in_curr_x = QLineEdit();    self.in_curr_x.setPlaceholderText("(m/s)")
        # self.in_curr_y = QLineEdit();    self.in_curr_y.setPlaceholderText("(m/s)")
        # self.in_curr_z = QLineEdit();    self.in_curr_z.setPlaceholderText("(m/s)")

        # _dv_env = QDoubleValidator(0.0, 1e12, 6)
        # for le in (self.in_dt, self.in_duration, self.in_curr_x, self.in_curr_y, self.in_curr_z):
        #     le.setValidator(_dv_env); le.setAlignment(Qt.AlignRight)

        # env_form.addRow("time step", self.in_dt)
        # env_form.addRow("total time", self.in_duration)
        # env_form.addRow("current x", self.in_curr_x)
        # env_form.addRow("current y", self.in_curr_y)
        # env_form.addRow("current z", self.in_curr_z)
        # form.addRow(self.grp_env)


        # 🔹 SpinBox → LineEdit + 콤마 표시
        _dv = QDoubleValidator(0.0, 1e12, 6)

        self.in_E = QLineEdit(w)
        self.in_rho = QLineEdit(w)
        self.in_d = QLineEdit(w)
        for le in (self.in_E, self.in_rho,):
            le.setValidator(_dv)
            le.setAlignment(Qt.AlignRight)
            le.setPlaceholderText("0")

            # 입력 시 콤마 자동 추가
            le.textChanged.connect(lambda text, line=le: self._format_with_comma(line, text))

        self.in_d.setValidator(_dv)
        self.in_d.setAlignment(Qt.AlignRight)
        self.in_d.setPlaceholderText("0")
        self.in_d.textChanged.connect(lambda text, line=self.in_d: self._format_with_comma(line, text))

        form.addRow("Young's modulus (Pa)", self.in_E)
        form.addRow("Density (kg/m³)", self.in_rho)
        form.addRow("Diameter (m)", self.in_d)

        # =========================
        #   NET 섹션 (박스)
        # =========================
        self.grp_net = QGroupBox("NET")
        self.grp_net.setCheckable(False)
        self.grp_net.setStyleSheet("""
        QGroupBox {
            font-weight: 600;
            border: 1px solid #bbb;
            border-radius: 6px;
            margin-top: 8px;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 4px;
        }
        """)
        net_form = QFormLayout(self.grp_net)

        self.in_mesh_length = QLineEdit()
        self.in_mesh_length.setPlaceholderText("(m)")
        self.in_bottom_net_weight = QLineEdit()
        self.in_bottom_net_weight.setPlaceholderText("(kg)")

        _dv = QDoubleValidator(0.0, 1e9, 4)
        self.in_mesh_length.setValidator(_dv)
        self.in_bottom_net_weight.setValidator(_dv)

        net_form.addRow("Mesh length", self.in_mesh_length)
        self.chk_bottom_net_center = QCheckBox()
        self.chk_bottom_net_center.setChecked(False)
        net_form.addRow("Bottom Net Center", self.chk_bottom_net_center)
        self.in_bottom_net_weight.setEnabled(False)
        self.chk_bottom_net_center.toggled.connect(
            lambda checked: self.in_bottom_net_weight.setEnabled(checked)
        )
        net_form.addRow("weight", self.in_bottom_net_weight)
        self.grp_net.setEnabled(False)
        form.addRow(self.grp_net)

        # =========================
        #   FLOATER 섹션 (박스)
        # =========================
        self.grp_floater = QGroupBox("FLOATER")
        self.grp_floater.setCheckable(False)
        self.grp_floater.setStyleSheet("""
        QGroupBox {
            font-weight: 600;
            border: 1px solid #bbb;
            border-radius: 6px;
            margin-top: 8px;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 4px;
        }
        """)
        flo_form = QFormLayout(self.grp_floater)

        self.in_floater_diam = QLineEdit()
        self.in_floater_diam.setPlaceholderText("(m)")
        self.in_floater_thk = QLineEdit()
        self.in_floater_thk.setPlaceholderText("(m)")
        self.in_floater_diam.setValidator(_dv)
        self.in_floater_thk.setValidator(_dv)

        flo_form.addRow("Floater diameter", self.in_floater_diam)
        flo_form.addRow("Floater thickness", self.in_floater_thk)
        self.grp_floater.setEnabled(False)
        form.addRow(self.grp_floater)

        # =========================
        #   BOTTOM RING 섹션 (박스)
        # =========================
        self.grp_bottom = QGroupBox("SINKER") 
        self.grp_bottom.setCheckable(False)
        self.grp_bottom.setStyleSheet("""
        QGroupBox {
            font-weight: 600;
            border: 1px solid #bbb;
            border-radius: 6px;
            margin-top: 8px;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 4px;
        }
        """)
        sk_form = QFormLayout(self.grp_bottom)

        # ✅ single sinker weight 입력 (kg)
        self.in_single_sinker_kg = QLineEdit()
        self.in_single_sinker_kg.setPlaceholderText("(kg)")
        self.in_single_sinker_kg.setValidator(_dv)
        self.in_single_sinker_kg.setAlignment(Qt.AlignRight)

        sk_form.addRow("single sinker weight", self.in_single_sinker_kg)

        # 기본은 비활성
        self.grp_bottom.setEnabled(False)
        form.addRow(self.grp_bottom)

        # =========================
        #   BUOYLINE 섹션 (박스)
        # =========================
        self.grp_buoy = QGroupBox("BUOYLINE")
        self.grp_buoy.setCheckable(False)
        self.grp_buoy.setStyleSheet("""
        QGroupBox {
            font-weight: 600;
            border: 1px solid #bbb;
            border-radius: 6px;
            margin-top: 8px;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 4px;
        }
        """)
        buoy_form = QFormLayout(self.grp_buoy)

        self.in_buoy_force = QLineEdit()
        self.in_buoy_force.setPlaceholderText("(N) per buoy node")
        self.in_buoy_force.setValidator(_dv)
        self.in_buoy_force.setAlignment(Qt.AlignRight)

        buoy_form.addRow("Buoyancy", self.in_buoy_force)

        # 기본은 비활성 (buoyline1/2 선택 시만 활성)
        self.grp_buoy.setEnabled(False)
        form.addRow(self.grp_buoy)

        btn_save = QPushButton("Assign & Save", w)

        def _save():
            # 현재 보이는 그룹 값 먼저 버퍼에 반영
            key_now = None
            if hasattr(self, "group_checks"):
                for k, cb in self.group_checks.items():
                    if cb.isChecked():
                        key_now = k
                        break
            if key_now:
                self._save_form_to_state(key_now)

            # 버퍼 내용 합쳐서 JSON 스키마로 재배치
            data = {
                "Floater": {}, "Net": {}, "Bottom_ring": {}, "Sinker": {},
                "Side_ropes": {}, "Bridle": {}, "Buoyline1": {},
                "Anchor_A": {}, "Anchor_B": {}, "Mooring_line": {}, "Mooring_frame": {}
            }
            
            # 🔸 ENVIRONMENT 값도 materials.json에 포함
            def _getf(le, default=None):
                try:
                    t = (le.text() if le else "").strip().replace(",", "")
                    return float(t) if t != "" else default
                except Exception:
                    return default

            dt_val       = _getf(self.in_dt)
            duration_val = _getf(self.in_duration)
            curx         = _getf(self.in_curr_x)
            cury         = _getf(self.in_curr_y)
            curz         = _getf(self.in_curr_z)

            data["simulation"] = {
                "dt": dt_val if dt_val is not None else 0.0,
                "duration": duration_val if duration_val is not None else 0.0,
                "current": [curx or 0.0, cury or 0.0, curz or 0.0]
            }

            if hasattr(self, "grp_wave"):
                data["simulation"]["wave"] = {
                    "enabled": bool(self.grp_wave.isChecked()),
                    "type": "regular",
                    "height": _getf(self.in_wave_height, 0.0) or 0.0,
                    "period": _getf(self.in_wave_period, 1.0) or 1.0,
                    "wavelength": _getf(self.in_wave_length, 10.0) or 10.0,
                    "depth": _getf(self.in_wave_depth, 1.0) or 1.0,
                    "direction_deg": _getf(self.in_wave_dir, 0.0) or 0.0,
                    "phase_deg": _getf(self.in_wave_phase, 0.0) or 0.0,
                    "z0": _getf(self.in_wave_z0, 0.0) or 0.0,
                }
        
            st = self._form_state

            if "floater" in st:
                f = st["floater"]
                # 값이 모두 존재할 때만 기록(부분 입력은 무시하고 싶지 않다면 그대로 기록)
                data["Floater"] = {
                    "FloaterYoungModule": f.get("FloaterYoungModule"),
                    "RHO": f.get("RHO"),
                    "Floaterdiameter": f.get("Floaterdiameter"),
                    "Floaterthickness": f.get("Floaterthickness"),
                }
            if "net twine" in st:
                n = st["net twine"]

                # 🔹 여기서 Sn 자동 계산
                D = n.get("twine_dia")
                L = n.get("Mesh_length")
                Sn = None
                if D is not None and L not in (None, 0):
                    ratio = D / L
                    Sn = 2 * ratio + 0.5 * (ratio ** 2)

                data["Net"] = {
                    "NetYoungModule": n.get("NetYoungModule"),
                    "RHO": n.get("RHO"),
                    "twine_dia": n.get("twine_dia"),
                    "Mesh_length": n.get("Mesh_length"),
                    "bottom_net_center": bool(self.chk_bottom_net_center.isChecked()),
                    "BottomNet_cen_weight": (
                        n.get("BottomNet_cen_weight") * 9.8
                        if self.chk_bottom_net_center.isChecked() and n.get("BottomNet_cen_weight") is not None
                        else 0.0
                    ),
                    "Sn": Sn,   # 🔥 여기 추가: D, L로부터 계산된 Sn
                }
            if "bottom ring" in st:
                b = st["bottom ring"]
                data["Bottom_ring"] = {
                    "BrYoungModule": b.get("BrYoungModule"),
                    "RHO": b.get("RHO"),
                    "Br_diameter": b.get("Br_diameter"),
                }
            if "bracket" in st:  # ✅ 추가
                b = st["bracket"]
                data["Bracket"] = {
                    "BrktYoungModule": b.get("BrktYoungModule"),
                    "RHO":             b.get("RHO"),
                    "Brkt_diameter":   b.get("Brkt_diameter"),
                }

            if hasattr(self, "group_checks") and "sinker" in self.group_checks:
                if self.group_checks["sinker"].isEnabled():  # sinker 그룹이 존재할 때만
                    try:
                        kg = float((self.in_single_sinker_kg.text() or "0").strip().replace(",", ""))
                    except Exception:
                        kg = 0.0

                    data["Sinker"] = {
                        "single_sinker_weight": kg * 9.8
                    }

            if "side ropes" in st:
                s = st["side ropes"]
                data["Side_ropes"] = {
                    "SideRopeYoungModule": s.get("SideRopeYoungModule"),
                    "RHO": s.get("RHO"),
                    "SideRope_diameter": s.get("SideRope_diameter"),
                }

            if "bridle line" in st:
                s = st["bridle line"]
                data["Bridle"] = {
                    "BridleYoungModule": s.get("BridleYoungModule"),
                    "RHO": s.get("RHO"),
                    "Bridle_diameter": s.get("Bridle_diameter"),
                }

            if "buoyline1" in st:
                s = st["buoyline1"]
                data["Buoyline1"] = {
                    "Buoy1YoungModule": s.get("Buoy1YoungModule"),
                    "RHO": s.get("RHO"),
                    "Buoy1_diameter": s.get("Buoy1_diameter"),
                    "BuoyForce": s.get("BuoyForce"),
                }

            # ✅ Buoyline2는 "있을 때만" 섹션 자체를 추가 저장
            if "buoyline2" in st:
                s = st["buoyline2"]
                data["Buoyline2"] = {
                    "Buoy2YoungModule": s.get("Buoy2YoungModule"),
                    "RHO": s.get("RHO"),
                    "Buoy2_diameter": s.get("Buoy2_diameter"),
                    "BuoyForce": s.get("BuoyForce"),
                }

            if "anchor line 1" in st:
                s = st["anchor line 1"]
                data["Anchor_A"] = {
                    "AnchorAYoungModule": s.get("AnchorAYoungModule"),
                    "RHO": s.get("RHO"),
                    "AnchorA_diameter": s.get("AnchorA_diameter"),
                }

            if "anchor line 2" in st:
                s = st["anchor line 2"]
                data["Anchor_B"] = {
                    "AnchorBYoungModule": s.get("AnchorBYoungModule"),
                    "RHO": s.get("RHO"),
                    "AnchorB_diameter": s.get("AnchorB_diameter"),
                }

            if "mooring line" in st:
                s = st["mooring line"]
                data["Mooring_line"] = {
                    "MooringYoungModule": s.get("MooringYoungModule"),
                    "RHO": s.get("RHO"),
                    "Mooring_diameter": s.get("Mooring_diameter"),
                }

            if "mooring frame" in st:
                s = st["mooring frame"]
                data["Mooring_frame"] = {
                    "MoorFrameYoungModule": s.get("MoorFrameYoungModule"),
                    "RHO": s.get("RHO"),
                    "MoorFrame_diameter": s.get("MoorFrame_diameter"),
                }

            # 🔸 Save 자동 경로 지정: v2024/Fish_Cage/asterinput/materials.py
            base = os.path.dirname(os.path.abspath(__file__))             # .../Fish_Cage/Fish_cage_simulation
            root = os.path.abspath(os.path.join(base, os.pardir))         # .../Fish_Cage
            aster = os.path.join(root, "asterinput")                      # .../Fish_Cage/asterinput
            os.makedirs(aster, exist_ok=True)

            save_path = os.path.join(aster, "materials.py")

            try:
                with open(save_path, "w", encoding="utf-8") as f:
                    f.write("materials = ")
                    f.write(repr(data))   # ✅ 파이썬 dict 형태로 저장 (True/False, None 등)
                    f.write("\n")

                QMessageBox.information(
                    self,
                    "Saved",
                    f"모든 그룹 재료 정보를 자동 저장했습니다.\n→ {save_path}"
                )

            except Exception as e:
                QMessageBox.critical(self, "저장 실패", f"파일 저장 중 오류: {e}")



        btn_save.clicked.connect(_save)
        form.addRow("", btn_save)

        parent_layout.addWidget(w)

    # 🔸 콤마 포맷 함수 추가 (클래스 내부 어디든)
    def _format_with_comma(self, line_edit, text):
        if not text:
            return
        # 소수 입력 진행 중이면(끝이 '.' 또는 '0'이고, 소수점 포함) 포맷하지 않음
        if '.' in text and (text.endswith('.') or text.endswith('0')):
            return

        clean = text.replace(",", "")
        if clean in ('-',):
            return
        try:
            num = float(clean)
        except ValueError:
            return

        if '.' in clean:
            fmt = f"{num:,.6f}".rstrip("0").rstrip(".")
        else:
            fmt = f"{int(num):,}"

        line_edit.blockSignals(True)
        line_edit.setText(fmt)
        line_edit.blockSignals(False)


    def _set_net_section_enabled(self, on: bool):
        if hasattr(self, "grp_net"):
            self.grp_net.setEnabled(bool(on))

    def _set_floater_section_enabled(self, on: bool):
        if hasattr(self, "grp_floater"):
            self.grp_floater.setEnabled(bool(on))

    def attach_environment_form(self, parent_layout):
        w = QWidget(self)
        root = QVBoxLayout(w)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        _dv_env = QDoubleValidator(-1e12, 1e12, 6)

        def _style_groupbox(gb):
            gb.setCheckable(False)
            gb.setStyleSheet("""
            QGroupBox {
                font-weight: 600;
                border: 1px solid #bbb;
                border-radius: 6px;
                margin-top: 8px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 10px;
                padding: 0 4px;
            }
            """)

        # =========================
        # TIME
        # =========================
        self.grp_time = QGroupBox("TIME")
        _style_groupbox(self.grp_time)
        time_form = QFormLayout(self.grp_time)

        self.in_dt = QLineEdit()
        self.in_dt.setPlaceholderText("(s)")
        self.in_duration = QLineEdit()
        self.in_duration.setPlaceholderText("(s)")

        for le in (self.in_dt, self.in_duration):
            le.setValidator(_dv_env)
            le.setAlignment(Qt.AlignRight)

        time_form.addRow("time step", self.in_dt)
        time_form.addRow("total time", self.in_duration)

        root.addWidget(self.grp_time)

        # =========================
        # CURRENT
        # =========================
        self.grp_current = QGroupBox("CURRENT")
        _style_groupbox(self.grp_current)
        cur_form = QFormLayout(self.grp_current)

        self.in_curr_x = QLineEdit()
        self.in_curr_x.setPlaceholderText("(m/s)")
        self.in_curr_y = QLineEdit()
        self.in_curr_y.setPlaceholderText("(m/s)")
        self.in_curr_z = QLineEdit()
        self.in_curr_z.setPlaceholderText("(m/s)")

        for le in (self.in_curr_x, self.in_curr_y, self.in_curr_z):
            le.setValidator(_dv_env)
            le.setAlignment(Qt.AlignRight)

        cur_form.addRow("current x", self.in_curr_x)
        cur_form.addRow("current y", self.in_curr_y)
        cur_form.addRow("current z", self.in_curr_z)

        root.addWidget(self.grp_current)

        # =========================
        # WAVE
        # =========================
        self.grp_wave = QGroupBox("WAVE")
        self.grp_wave.setCheckable(True)
        self.grp_wave.setChecked(False)

        wave_form = QFormLayout(self.grp_wave)

        self.in_wave_height = QLineEdit()
        self.in_wave_height.setPlaceholderText("(m)")

        self.in_wave_period = QLineEdit()
        self.in_wave_period.setPlaceholderText("(s)")

        self.in_wave_length = QLineEdit()
        self.in_wave_length.setPlaceholderText("(m)")

        self.in_wave_depth = QLineEdit()
        self.in_wave_depth.setPlaceholderText("(m)")

        self.in_wave_dir = QLineEdit()
        self.in_wave_dir.setPlaceholderText("(deg)")

        self.in_wave_phase = QLineEdit()
        self.in_wave_phase.setPlaceholderText("(deg)")

        self.in_wave_z0 = QLineEdit()
        self.in_wave_z0.setPlaceholderText("(m)")

        _dv_wave = QDoubleValidator(-1e12, 1e12, 6)

        for le in (
            self.in_wave_height,
            self.in_wave_period,
            self.in_wave_length,
            self.in_wave_depth,
            self.in_wave_dir,
            self.in_wave_phase,
            self.in_wave_z0,
        ):
            le.setValidator(_dv_wave)
            le.setAlignment(Qt.AlignRight)

        wave_form.addRow("height", self.in_wave_height)
        wave_form.addRow("period", self.in_wave_period)
        wave_form.addRow("wavelength", self.in_wave_length)
        wave_form.addRow("depth", self.in_wave_depth)
        wave_form.addRow("direction", self.in_wave_dir)
        wave_form.addRow("phase", self.in_wave_phase)
        wave_form.addRow("z0", self.in_wave_z0)

        root.addWidget(self.grp_wave)
        root.addStretch(1)

        parent_layout.addWidget(w)

def main():
    app = QApplication(sys.argv)
    window = MainWindow()   # ← 존재하는 클래스명으로 호출
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
