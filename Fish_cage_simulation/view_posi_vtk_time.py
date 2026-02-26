import sys
import os
import glob
import re
import numpy as np

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QFileDialog, QToolBar, QAction, QMessageBox, QLabel, QSlider
)
from PyQt5.QtCore import Qt
from PyQt5 import QtCore

from vtk.qt.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
import vtk


def extract_time_from_name(path: str) -> float:
    """
    파일 이름에서 숫자(시간)를 추출해서 정렬에 사용.
    예: posi_0.0.txt, posi_9.9.txt, posi_10.2.txt → 0.0, 9.9, 10.2
    숫자가 없으면 0.0 리턴.
    """
    name = os.path.basename(path)
    m = re.search(r"([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", name)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return 0.0
    return 0.0


class PosiVTKViewer(QMainWindow):
    def __init__(self, path=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("posi VTK 3D Viewer (time series)")
        self.resize(1000, 800)

        # --- 상태 변수 (시간 시퀀스 관리) ---
        self.files = []         # 시간 순서대로 정렬된 posi 파일 리스트
        self.current_index = -1 # 현재 보고 있는 파일 인덱스
        self.posi_dir = None

        # --- meshinfo 관련 상태 ---
        self.meshinfo = None
        self.meshinfo_path = None
        self.lines_netting = None
        self.lines_pipe_top = None
        self.lines_bracket = None
        self.lines_bottom_ring = None
        self.node_actor = None
        self.nodes_visible = True
        self.lines_side_rope = None
        self.lines_bridle    = None
        self.lines_buoy      = None
        self.lines_anchor_a  = None
        self.lines_anchor_b  = None
        self.lines_moorfrm   = None
        self._camera_initialized = False

        # --- VTK 위젯 생성 ---
        self.vtk_widget = QVTKRenderWindowInteractor(self)
        self.setCentralWidget(self.vtk_widget)

        self.renderer = vtk.vtkRenderer()
        self.vtk_widget.GetRenderWindow().AddRenderer(self.renderer)

        # 인터랙터 스타일 (마우스 회전/줌)
        self.interactor = self.vtk_widget.GetRenderWindow().GetInteractor()
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)

        # Axes 표시
        self.axes = vtk.vtkAxesActor()
        self.axes_widget = vtk.vtkOrientationMarkerWidget()
        self.axes_widget.SetOrientationMarker(self.axes)
        self.axes_widget.SetInteractor(self.interactor)
        self.axes_widget.SetEnabled(1)
        self.axes_widget.SetOutlineColor(0.3, 0.3, 0.3)
        self.axes_widget.SetViewport(0.0, 0.0, 0.2, 0.2)

        # 툴바
        toolbar = QToolBar("MainToolbar")
        self.addToolBar(toolbar)

        act_open = QAction("Open posi", self)
        act_open.triggered.connect(self.open_file_dialog)
        toolbar.addAction(act_open)

        act_prev = QAction("Prev", self)
        act_prev.triggered.connect(self.show_prev)
        toolbar.addAction(act_prev)

        act_next = QAction("Next", self)
        act_next.triggered.connect(self.show_next)
        toolbar.addAction(act_next)

        act_toggle_nodes = QAction("Nodes ON/OFF", self)
        act_toggle_nodes.triggered.connect(self.toggle_nodes)
        toolbar.addAction(act_toggle_nodes)
        # --- 자동재생용 타이머 ---
        self.timer = QtCore.QTimer()
        self.timer.setInterval(0)  # 150ms = 약 6~7fps

        act_start = QAction("Start", self)
        act_start.triggered.connect(self.start_auto)
        toolbar.addAction(act_start)

        act_stop = QAction("Stop", self)
        act_stop.triggered.connect(self.stop_auto)
        toolbar.addAction(act_stop)

        # 상태표시 (현재 시간/인덱스)
        self.status_label = QLabel("")
        self.statusBar().addPermanentWidget(self.status_label)

        # --- Time slider (0 ~ last index) ---
        self.time_slider = QSlider(Qt.Horizontal)
        self.time_slider.setMinimum(0)
        self.time_slider.setMaximum(0)
        self.time_slider.setValue(0)
        self.time_slider.setFixedWidth(260)  # 취향: 길이 조절
        self.time_slider.valueChanged.connect(self.on_slider_changed)
        self.statusBar().addPermanentWidget(self.time_slider)

        # 초기 파일 있으면 바로 로드
        if path is not None and os.path.exists(path):
            self.setup_sequence_and_show(path)

        self.interactor.Initialize()

    # ----------------- meshinfo.py 찾기 & 로딩 -----------------
    def _find_meshinfo_from_posi(self, posi_path: str) -> str | None:
        """
        posi 파일이 있는 위치에서 상위 폴더들을 타고 올라가면서
        meshinfo.py 를 찾는다.

        보통 구조:
        v2024/Fish_Cage/asterinput/pythonOutput/posi_*.txt
        → v2024/Fish_Cage/asterinput/meshinfo.py
        """
        cur = os.path.abspath(os.path.dirname(posi_path))
        while True:
            # 1) 현재 폴더에 meshinfo.py 가 있는지 우선 확인
            candidate_here = os.path.join(cur, "meshinfo.py")
            if os.path.isfile(candidate_here):
                return candidate_here

            # 2) asterinput 서브 폴더를 직접 찾는 경우
            candidate_aster = os.path.join(cur, "asterinput", "meshinfo.py")
            if os.path.isfile(candidate_aster):
                return candidate_aster

            parent = os.path.dirname(cur)
            if parent == cur:
                break
            cur = parent
        return None
    
    def _sync_slider(self):
        """현재 index/파일수와 슬라이더 상태 동기화"""
        if not hasattr(self, "time_slider"):
            return
        n = len(self.files) if self.files else 1
        maxi = max(0, n - 1)
        cur = 0 if self.current_index < 0 else min(self.current_index, maxi)

        self.time_slider.blockSignals(True)
        self.time_slider.setMinimum(0)
        self.time_slider.setMaximum(maxi)
        self.time_slider.setValue(cur)
        self.time_slider.blockSignals(False)

    def on_slider_changed(self, value: int):
        if not self.files:
            return
        self.refresh_sequence()
        self.show_by_index(int(value))

    def _ensure_meshinfo_loaded(self, posi_path: str):
        """
        한 번만 meshinfo.py를 읽어서
        Lines_netting / Lines_pipe_top / Line_braket 정보를 준비.
        """
        if self.meshinfo is not None:
            return

        mpath = self._find_meshinfo_from_posi(posi_path)
        if not mpath:
            print("※ meshinfo.py를 찾지 못했습니다. (net/floater/bracket 선은 표시되지 않습니다.)")
            self.meshinfo = {}
            self.lines_netting = None
            self.lines_pipe_top = None
            self.lines_bracket = None
            return

        scope = {}
        try:
            with open(mpath, "r", encoding="utf-8") as f:
                code = f.read()
            exec(code, scope)
        except Exception as e:
            print(f"※ meshinfo.py 로딩 실패: {e}")
            self.meshinfo = {}
            self.lines_netting = None
            self.lines_pipe_top = None
            self.lines_bracket = None
            return

        meshinfo = scope.get("meshinfo", None)
        if not isinstance(meshinfo, dict):
            print("※ meshinfo.py 안에 meshinfo 딕셔너리가 없습니다.")
            self.meshinfo = {}
            self.lines_netting = None
            self.lines_pipe_top = None
            self.lines_bracket = None
            return

        self.meshinfo = meshinfo
        self.meshinfo_path = mpath

        # Lines_* 는 MED 노드 번호(1-베이스)를 쓰므로
        # VTK points 인덱스에 쓸 때는 -1 해서 0-베이스로 쓴다.
        self.lines_netting = meshinfo.get("Lines_netting") or []
        self.lines_pipe_top = meshinfo.get("Lines_pipe_top") or []
        self.lines_bracket = meshinfo.get("Line_braket") or []
        self.lines_bottom_ring = meshinfo.get("Lines_bottom_ring") or []

        self.lines_side_rope = meshinfo.get("Lines_side_rope") or []
        self.lines_bridle    = meshinfo.get("Lines_bridle") or []
        self.lines_buoy      = meshinfo.get("Lines_buoy") or []
        self.lines_anchor_a  = meshinfo.get("Lines_anchor_A") or []
        self.lines_anchor_b  = meshinfo.get("Lines_anchor_B") or []
        self.lines_moorfrm   = meshinfo.get("Lines_mooring_frame") or []

        print(f"[meshinfo] {mpath} 로드 완료")
        print(f"  Lines_netting : {len(self.lines_netting)}")
        print(f"  Lines_pipe_top: {len(self.lines_pipe_top)}")
        print(f"  Line_braket   : {len(self.lines_bracket)}")

    # ----------------- 파일 열기 & 시퀀스 구성 -----------------
    def open_file_dialog(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open posi text file",
            "",
            "Text Files (*.txt *.dat);;All Files (*)"
        )
        if not path:
            return
        self.setup_sequence_and_show(path)

    def refresh_sequence(self):
        """posi_dir에서 새로 생긴 posi*.txt를 감지해 self.files를 갱신"""
        if not self.posi_dir or not os.path.isdir(self.posi_dir):
            return False

        candidates = glob.glob(os.path.join(self.posi_dir, "posi*.txt"))
        if not candidates:
            return False

        candidates.sort(key=extract_time_from_name)

        # 변화 없으면 종료
        if len(candidates) == len(self.files):
            return False

        # 현재 보고 있는 파일 시간 유지하면서 index 재설정
        cur_t = None
        if 0 <= self.current_index < len(self.files):
            cur_t = extract_time_from_name(self.files[self.current_index])

        self.files = candidates

        if cur_t is not None:
            times = [extract_time_from_name(f) for f in self.files]
            self.current_index = min(range(len(times)), key=lambda i: abs(times[i] - cur_t))
        else:
            self.current_index = min(max(self.current_index, 0), len(self.files) - 1)

        self._sync_slider()
        return True

    def setup_sequence_and_show(self, path: str):
        """
        1) 선택한 path를 기준으로 같은 폴더의 posi*.txt 파일들을 전부 찾고
        2) 시간(파일명 속 숫자) 기준으로 정렬한 뒤
        3) 그 리스트에서 path 위치를 current_index로 설정하고 화면에 표시
        """
        directory = os.path.dirname(path)
        self.posi_dir = directory

        # 같은 폴더 내 posi*.txt 파일 모두 찾기
        candidates = glob.glob(os.path.join(directory, "posi*.txt"))
        if not candidates:
            # 그래도 현재 파일만 보여주기
            self.files = [path]
            self.current_index = 0
            self.load_posi(path)
            return

        # 시간 값으로 정렬
        candidates.sort(key=extract_time_from_name)

        self.files = candidates
        # 현재 파일이 리스트에서 몇 번째인지 찾기
        try:
            self.current_index = self.files.index(path)
        except ValueError:
            # 혹시 패턴이 안 맞아서 못 찾으면, 가장 가까운 시간으로 위치 설정
            t_target = extract_time_from_name(path)
            times = [extract_time_from_name(f) for f in self.files]
            idx = min(range(len(times)), key=lambda i: abs(times[i] - t_target))
            self.current_index = idx

        self.show_by_index(self.current_index)
        self._sync_slider()

    def toggle_nodes(self):
        if self.actor_points is None:
            return

        self.nodes_visible = not self.nodes_visible
        self.actor_points.SetVisibility(1 if self.nodes_visible else 0)
        self.vtk_widget.GetRenderWindow().Render()

    # ----------------- 시퀀스 이동 -----------------
    def show_prev(self):
        self.refresh_sequence()
        if not self.files:
            return

        if self.current_index <= 0:
            self.current_index = 0
            self.show_by_index(self.current_index)
            return

        self.current_index -= 1
        self.show_by_index(self.current_index)

    def show_next(self):
        if not self.files:
            return

        old_last = len(self.files) - 1
        old_idx = self.current_index

        self.refresh_sequence()  # 새 파일 감지
        new_last = len(self.files) - 1

        # (1) 예전엔 마지막이었는데, 새 파일이 생겼으면 다음으로 이동
        if old_idx == old_last and new_last > old_last:
            self.current_index = min(old_idx + 1, new_last)
            self.show_by_index(self.current_index)
            return

        # (2) 새 파일이 없으면 "끝에서 멈춤" (wrap-around 금지)
        if self.current_index >= new_last:
            self.current_index = new_last
            self.show_by_index(self.current_index)

            # 자동재생 중이면 여기서 멈추고 싶으면 주석 해제
            # if self.timer.isActive():
            #     self.stop_auto()
            return

        # (3) 일반적으로 다음 인덱스
        self.current_index += 1
        self.show_by_index(self.current_index)

    def show_by_index(self, index: int):
        if not self.files:
            return
        index = max(0, min(index, len(self.files) - 1))
        path = self.files[index]
        self.current_index = index
        self.load_posi(path)
        self._sync_slider()
    # ----------------- 실제 posi 파일 로딩 & VTK 표시 -----------------
    def load_posi(self, path):
        try:
            arr = np.loadtxt(path)
        except Exception as e:
            QMessageBox.critical(self, "Load error", f"파일을 읽을 수 없습니다.\n{e}")
            return

        arr = np.asarray(arr)

        # (N,3) 강제
        if arr.ndim == 1:
            if arr.size % 3 != 0:
                QMessageBox.critical(self, "Shape error", "데이터를 (N,3)로 reshape 할 수 없습니다.")
                return
            arr = arr.reshape(-1, 3)

        if arr.shape[1] != 3:
            QMessageBox.critical(self, "Shape error", f"(N,3) 형식이 아닙니다: {arr.shape}")
            return

        # meshinfo(선 정보) 준비
        self._ensure_meshinfo_loaded(path)

        x, y, z = arr[:, 0], arr[:, 1], arr[:, 2]
        n_nodes = arr.shape[0]

        # VTK Points 생성
        points = vtk.vtkPoints()
        for xi, yi, zi in zip(x, y, z):
            points.InsertNextPoint(float(xi), float(yi), float(zi))

        # --- 점 표시용 ---
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)

        glyph_filter = vtk.vtkVertexGlyphFilter()
        glyph_filter.SetInputData(polydata)
        glyph_filter.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(glyph_filter.GetOutputPort())

        self.actor_points = vtk.vtkActor()
        self.actor_points.SetMapper(mapper)
        self.actor_points.GetProperty().SetPointSize(4)


        self.actor_points.SetVisibility(1 if self.nodes_visible else 0)
        # ==========================================================
        # 🔥 Line element + 노드 인장력(N) → line 평균 인장력 컬러
        # ==========================================================
        
        cam_state = self._get_camera_state() if self._camera_initialized else None

        # renderer 초기화
        self.renderer.RemoveAllViewProps()
        self.renderer.AddActor(self.actor_points)

        # ------------------------------
        # 1) 노드 인장력 N 로드
        # ------------------------------
        base = os.path.dirname(os.path.dirname(path))
        N_dir = os.path.join(base, "N")

        N = None
        if os.path.isdir(N_dir):
            files = glob.glob(os.path.join(N_dir, "N*.txt"))
            if files:
                files.sort(key=lambda p: abs(extract_time_from_name(p) - extract_time_from_name(path)))
                try:
                    N = np.loadtxt(files[0], dtype=float)
                except Exception as e:
                    print(f"[WARN] N load failed: {e}")

        if N is None or len(N) != n_nodes:
            print(f"[WARN] N invalid: len(N)={None if N is None else len(N)}, n_nodes={n_nodes}")
            N = None

        # ------------------------------
        # 2) line cell 생성 (net + pipe + bracket)
        # ------------------------------
        cell_array = vtk.vtkCellArray()
        line_pairs = []

        for group in (
            self.lines_netting,
            self.lines_pipe_top,
            self.lines_bracket,
            self.lines_bottom_ring,
            self.lines_side_rope,
            self.lines_bridle,
            self.lines_buoy,
            self.lines_anchor_a,
            self.lines_anchor_b,
            self.lines_moorfrm,
        ):
            if not group:
                continue
            for a_med, b_med in group:
                i = int(a_med) - 1
                j = int(b_med) - 1
                if 0 <= i < n_nodes and 0 <= j < n_nodes:
                    line = vtk.vtkLine()
                    line.GetPointIds().SetId(0, i)
                    line.GetPointIds().SetId(1, j)
                    cell_array.InsertNextCell(line)
                    line_pairs.append((i, j))

        # ------------------------------
        # 3) PolyData
        # ------------------------------
        line_poly = vtk.vtkPolyData()
        line_poly.SetPoints(points)
        line_poly.SetLines(cell_array)

        # ------------------------------
        # 4) line 인장력 = (노드 i + 노드 j)/2
        # ------------------------------
        if N is not None:
            ncell = line_poly.GetNumberOfCells()

            arr = vtk.vtkFloatArray()
            arr.SetName("Tension")
            arr.SetNumberOfTuples(ncell)

            for k, (i, j) in enumerate(line_pairs):
                Tline = 0.5 * (float(N[i]) + float(N[j]))
                arr.SetValue(k, Tline)

            line_poly.GetCellData().SetScalars(arr)
            line_poly.GetCellData().SetActiveScalars("Tension")
            line_poly.Modified()

        # ------------------------------
        # 5) Mapper + 컬러맵
        # ------------------------------
        mapper_l = vtk.vtkPolyDataMapper()
        mapper_l.SetInputData(line_poly)

        if N is not None:
            mapper_l.SetScalarModeToUseCellData()
            mapper_l.SetColorModeToMapScalars()
            mapper_l.ScalarVisibilityOn()

            vmin, vmax = arr.GetRange()
            if abs(vmax - vmin) < 1e-12:
                vmax = vmin + 1e-12
            mapper_l.SetScalarRange(vmin, vmax)

            lut = vtk.vtkLookupTable()
            lut.SetNumberOfTableValues(256)
            lut.Build()

            for i in range(256):
                x = i / 255.0

                if x < 0.5:
                    # Blue → White
                    t = x / 0.5
                    r = 0.9 * t + 0.10
                    g = 0.9 * t + 0.10
                    b = 10.00
                else:
                    # White → Red
                    t = (x - 0.5) / 0.5
                    r = 1.00
                    g = 0.85 * (1 - t) + 0.15
                    b = 0.85 * (1 - t) + 0.15

                lut.SetTableValue(i, r, g, b, 1.0)

            mapper_l.SetLookupTable(lut)
        else:
            mapper_l.ScalarVisibilityOff()

        # ------------------------------
        # 6) Actor
        # ------------------------------
        actor_lines = vtk.vtkActor()
        actor_lines.SetMapper(mapper_l)
        actor_lines.GetProperty().SetLineWidth(2.0)

        self.renderer.AddActor(actor_lines)
        
        self.renderer.SetBackground(0.1, 0.1, 0.1)
        if not self._camera_initialized:
            self.renderer.ResetCamera()
            self._camera_initialized = True

        # 상태 표시 (시간 + 인덱스)
        t = extract_time_from_name(path)
        idx_text = f"{self.current_index + 1}/{len(self.files)}" if self.files else "1/1"
        self.status_label.setText(f"t = {t:.3f}  |  {idx_text}")

        self.setWindowTitle(f"posi VTK 3D Viewer - {os.path.basename(path)}")
        self.vtk_widget.GetRenderWindow().Render()

    # ----------------- 자동재생 기능 -----------------
    def start_auto(self):
        if not self.files:
            return
        try:
            self.timer.timeout.disconnect(self.show_next)
        except Exception:
            pass
        self.timer.timeout.connect(self.show_next)
        self.timer.start()

    def stop_auto(self):
        self.timer.stop()
        # 중복 연결 방지
        try:
            self.timer.timeout.disconnect(self.show_next)
        except Exception:
            pass

    def _get_camera_state(self):
        cam = self.renderer.GetActiveCamera()
        return (
            cam.GetPosition(),
            cam.GetFocalPoint(),
            cam.GetViewUp(),
            cam.GetParallelScale(),
            cam.GetViewAngle(),
            cam.GetClippingRange(),
            cam.GetParallelProjection(),
        )

    def _set_camera_state(self, st):
        if st is None:
            return
        cam = self.renderer.GetActiveCamera()
        pos, foc, up, pscale, vangle, clip, pproj = st
        cam.SetPosition(*pos)
        cam.SetFocalPoint(*foc)
        cam.SetViewUp(*up)
        cam.SetParallelScale(pscale)
        cam.SetViewAngle(vangle)
        cam.SetClippingRange(*clip)
        cam.SetParallelProjection(pproj)

def main():
    app = QApplication(sys.argv)

    path = sys.argv[1] if len(sys.argv) > 1 else None

    win = PosiVTKViewer(path=path)
    win.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    import sys, os, glob

    app = QApplication(sys.argv)

    # -----------------------------------------
    # 🔹 자동 posi 디렉토리 설정
    # -----------------------------------------
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    POSI_DIR = os.path.abspath(
        os.path.join(
            BASE_DIR,
            "..",              # Fish_cage
            "asterinput",
            "pythonOutput",
            "posi"
        )
    )

    viewer = None

    if os.path.isdir(POSI_DIR):
        files = glob.glob(os.path.join(POSI_DIR, "posi*.txt"))
        if files:
            # 시간 기준 정렬
            files.sort(key=extract_time_from_name)
            viewer = PosiVTKViewer(path=files[0])
        else:
            viewer = PosiVTKViewer()
    else:
        viewer = PosiVTKViewer()

    viewer.show()
    sys.exit(app.exec_())
