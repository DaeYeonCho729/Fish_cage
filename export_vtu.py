# export_vtu.py
# GUI 환경에서 Salome 없이 .vtu를 생성합니다.
# 로직은 export_med.py와 동일한 섹션/그룹 규칙을 따라 polydata를 구성합니다.

import os, json, glob
import math
import vtk

BASE_DIR  = "Fish_Cage"
LINKS_DIR = os.path.join(BASE_DIR, "links")

def jload(path, default=None):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {} if default is None else default

def export_all_to_vtu(outpath="Fish_Cage/saves/fish_cage.vtu"):
    """
    export_med.py와 같은 조립 순서를 따릅니다.
      1) net_points.json : 노드 추가, z==0 → floater node / 그 외 → net node
                           엣지: floater 노드끼리만 'floater', 나머지 'net twine'
      2) float_temp.json : 'floater' 엣지/노드
      3) bottom_collar_temp.json : baseN shift 반영, 'bottom ring' 엣지/노드
      4) mooring_temp.json : 'mooring frame node'
      5) links/*.json : baseN shift 반영, 그룹 매핑(사이드로프/브라이들/부이/앵커 등)
    결과는 PolyData(.vtu)로 저장, Edge(Cell)에는 문자열 배열 'edge_group'을 설정.
    Point에는 대분류 문자열 'node_group'을 설정(중복그룹은 대표 하나로 표기).
    """
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    # 누적 컨테이너
    points = []          # list[[x,y,z], ...]
    edges  = []          # list[(i,j), ...]  (글로벌 인덱스 기준)
    edge_group_name = [] # 각 edge에 대응되는 그룹 이름 (len == len(edges))

    # 포인트/엣지 그룹(이름→인덱스 모음)
    node_groups = { }    # name -> set(gi)
    # ▼ 겹침 제거용 캐시
    bottom_ring_nodes = set()
    bottom_ring_edges = set()   # (min(i,j), max(i,j))로 저장

    # ---------------- 1) Net ----------------
    net = jload(os.path.join(BASE_DIR, "net_points.json"), {})
    net_pts   = net.get("points", []) or []
    net_edges = net.get("edges", [])  or []

    baseN = len(points)
    # 포인트 추가 + 분류
    floater_nodes_from_net, net_nodes_from_net = set(), set()
    for i, p in enumerate(net_pts):
        x, y, z = float(p[0]), float(p[1]), float(p[2])
        points.append([x, y, z])
        gi = baseN + i
        if abs(z) < 1e-8:
            floater_nodes_from_net.add(gi)
        else:
            net_nodes_from_net.add(gi)

    # 엣지 추가(분류 규칙 동일)
    def _map_net_idx(idx: int) -> int:
        i = int(idx)
        return (baseN + i) if (0 <= i < len(net_pts)) else i

    for e in net_edges:
        i, j = int(e[0]), int(e[1])
        gi, gj = _map_net_idx(i), _map_net_idx(j)
        if not (0 <= gi < len(points) and 0 <= gj < len(points)):
            continue
        edges.append((gi, gj))
        if (gi in floater_nodes_from_net) and (gj in floater_nodes_from_net):
            edge_group_name.append("floater")
        else:
            edge_group_name.append("net twine")

    if floater_nodes_from_net:
        node_groups.setdefault("floater node", set()).update(floater_nodes_from_net)
    if net_nodes_from_net:
        node_groups.setdefault("net node", set()).update(net_nodes_from_net)

    # ---------------- 2) Floating collar ----------------
    flo = jload(os.path.join(BASE_DIR, "float_temp.json"), {})
    flo_pts   = flo.get("points", []) or []
    flo_edges = flo.get("edges", [])  or []
    if flo_pts or flo_edges:
        secN = len(points)
        for p in flo_pts:
            points.append([float(p[0]), float(p[1]), float(p[2])])
        if flo_pts:
            node_groups.setdefault("floater node", set()).update(range(secN, secN + len(flo_pts)))
        for (a, b) in flo_edges:
            edges.append((secN + int(a), secN + int(b)))
            edge_group_name.append("floater")

    # ---------------- 3) Bottom collar (baseN shift) ----------------
    bot = jload(os.path.join(BASE_DIR, "bottom_collar_temp.json"), {})
    bot_pts   = bot.get("points_new", []) or bot.get("points", []) or []
    bot_edges = bot.get("ring_edges", []) or bot.get("edges", []) or []
    bot_ring_indices = bot.get("ring_indices", []) or []
    bot_baseN = bot.get("baseN", None)

    if bot_pts or bot_edges or bot_ring_indices:
        curr_before = len(points)
        def _shift(n: int) -> int:
            n = int(n)
            if bot_baseN is None:
                return n
            shift = curr_before - int(bot_baseN)
            return n + shift if n >= int(bot_baseN) else n

        # 섹션 포인트(추가 새 포인트) 먼저 push
        extra_start = len(points)
        for p in bot_pts:
            points.append([float(p[0]), float(p[1]), float(p[2])])

        # 링 노드 그룹(글로벌 인덱스, 시프트 반영)
        if bot_ring_indices:
            ring_shifted = [_shift(i) for i in bot_ring_indices]
            bottom_ring_nodes.update(ring_shifted)
            node_groups.setdefault("bottom ring node", set()).update(ring_shifted)

        # 엣지 (저장 당시 글로벌 인덱스가 들어있을 수 있어 시프트)
        if bot_edges:
            for (a, b) in bot_edges:
                aa, bb = _shift(a), _shift(b)
                if 0 <= aa < len(points) and 0 <= bb < len(points):
                    edges.append((aa, bb))
                    bottom_ring_edges.add((min(aa, bb), max(aa, bb)))
                    edge_group_name.append("bottom ring")

    # ---------------- 4) Mooring frame ----------------
    moor = jload(os.path.join(BASE_DIR, "mooring_temp.json"), {})
    moor_pts   = moor.get("points", []) or []
    moor_edges = moor.get("edges", [])  or []
    if moor_pts or moor_edges:
        secN = len(points)
        for p in moor_pts:
            points.append([float(p[0]), float(p[1]), float(p[2])])
        node_groups.setdefault("mooring frame node", set()).update(range(secN, secN + len(moor_pts)))
        for (a, b) in moor_edges:
            edges.append((secN + int(a), secN + int(b)))
            edge_group_name.append("mooring frame")

    # ---------------- 5) Links/*.json (baseN shift) ----------------
    key2edge = {
        "SIDE_ROPE_MAIN": "side ropes",
        "SIDE_ROPE_SUB":  "side ropes",
        "BRIDLE_LINE":    "bridle line",
        "BUOY_TETHER":    "buoyline1",
        "BUOY_LINE":      "buoyline2",
        "ANCHOR_LINE_A":  "anchor line 1",
        "ANCHOR_LINE_B":  "anchor line 2",
        "DISTANCE_ROPE":  "mooring line",
        "FLOATER_BRACKET":"bracket",
    }

    if os.path.isdir(LINKS_DIR):
        for path in sorted(glob.glob(os.path.join(LINKS_DIR, "*.json"))):
            data = jload(path, {})
            extra_pts = data.get("extra_points", []) or []
            polylines = data.get("polylines", []) or []
            saved_baseN = data.get("baseN", None)

            curr_before = len(points)
            def _S(n: int) -> int:
                n = int(n)
                if saved_baseN is None:
                    return n
                shift = curr_before - int(saved_baseN)
                return n + shift if n >= int(saved_baseN) else n

            # 섹션 extra 포인트 추가 (로컬 0..N-1)
            secN = len(points)
            for p in extra_pts:
                points.append([float(p[0]), float(p[1]), float(p[2])])

            # 폴리라인 -> 연속 엣지 + 그룹
            for pl in polylines:
                raw = [int(x) for x in (pl.get("nodes", []) or [])]
                if len(raw) < 2:
                    continue
                ns = [_S(x) for x in raw]
                label = key2edge.get(pl.get("group", ""), None)
                for i in range(len(ns) - 1):
                    a, b = ns[i], ns[i + 1]
                    if 0 <= a < len(points) and 0 <= b < len(points):
                        edges.append((a, b))
                        edge_group_name.append(label or "links")

    # ==== De-duplicate: remove net twine edges that overlap bottom ring ====
    if bottom_ring_nodes or bottom_ring_edges:
        keep_edges = []
        keep_groups = []
        for (a, b), gname in zip(edges, edge_group_name):
            if gname == "net twine":
                pair = (min(a, b), max(a, b))
                # 조건1) bottom ring과 동일한 엣지, 또는
                # 조건2) 두 끝점이 모두 bottom ring 노드에 속하면 → 제거
                if pair in bottom_ring_edges or (a in bottom_ring_nodes and b in bottom_ring_nodes):
                    continue  # drop this edge
            keep_edges.append((a, b))
            keep_groups.append(gname)
        edges = keep_edges
        edge_group_name = keep_groups

    # --------- all node 그룹(편의) ---------
    node_groups.setdefault("all node", set()).update(range(len(points)))

    # ================= VTK PolyData 구성 =================
    vtk_pts = vtk.vtkPoints()
    for x, y, z in points:
        vtk_pts.InsertNextPoint(float(x), float(y), float(z))

    cells = vtk.vtkCellArray()
    for (a, b) in edges:
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, int(a))
        line.GetPointIds().SetId(1, int(b))
        cells.InsertNextCell(line)

    poly = vtk.vtkPolyData()
    poly.SetPoints(vtk_pts)
    poly.SetLines(cells)

    # --- Edge(Cell) Data: edge_group (string) ---
    arr_edge_group = vtk.vtkStringArray()
    arr_edge_group.SetName("edge_group")
    for name in edge_group_name:
        arr_edge_group.InsertNextValue(str(name) if name is not None else "")
    poly.GetCellData().AddArray(arr_edge_group)

    # --- Point Data: node_group (대표 하나) ---
    # 대표 우선순위: mooring frame node > bottom ring node > floater node > net node > all node
    priority = ["mooring frame node", "bottom ring node", "floater node", "net node", "all node"]
    label_of = [""] * len(points)
    for name in priority:
        ids = node_groups.get(name, set())
        for gi in ids:
            if 0 <= gi < len(points):
                if not label_of[gi]:
                    label_of[gi] = name
    arr_node_group = vtk.vtkStringArray()
    arr_node_group.SetName("node_group")
    for s in label_of:
        arr_node_group.InsertNextValue(s)
    poly.GetPointData().AddArray(arr_node_group)

    # 저장
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(outpath)
    writer.SetInputData(poly)
    writer.SetDataModeToBinary()
    if writer.Write() == 0:
        raise RuntimeError("VTU write failed")
    
        # ------------------------------
    # meshinfo.json 생성
    # ------------------------------
    try:
        # net_points.json에서 panel과 z 정보 로드
        net_json_path = os.path.join(BASE_DIR, "net_points.json")
        net_raw = {}
        if os.path.isfile(net_json_path):
            with open(net_json_path, "r", encoding="utf-8") as f:
                net_raw = json.load(f)
        z_float  = float(net_raw.get("z_float", 0.0))
        surfs_from_net = net_raw.get("surfs_netting", []) or []

        # 그룹 배열 이름이 코드 상에서 무엇이든, 엣지 루프를 만들 때 함께 만든
        #  - edges           : [[i,j], ...] (global index)
        #  - edge_group_name : ["floater"/"net twine"/"bottom ring"/...]
        # 가 이미 존재해야 한다. (export_vtu.py의 기존 로직이 생성)
        def _approx(a, b, tol=1e-6): return abs(a - b) <= tol

        # Lines_pipe_top: floater 그룹 & z≈z_float
        lines_pipe_top = []
        for (a, b), gname in zip(edges, edge_group_name):
            if gname == "floater":
                za = float(points[a][2]); zb = float(points[b][2])
                if _approx(za, z_float) and _approx(zb, z_float):
                    lines_pipe_top.append([int(a), int(b)])

        # Lines_netting: 'net twine' 그룹 전체
        lines_netting = [[int(a), int(b)]
                         for (a, b), gname in zip(edges, edge_group_name)
                         if gname == "net twine"]

        # Nodes: 전체 좌표
        nodes_xyz = [[float(x), float(y), float(z)] for (x, y, z) in points]

        # numberOf*
        n_pipe  = len(lines_pipe_top)
        n_twine = len(lines_netting)
        n_surfs = len(surfs_from_net)

        # weight_NT 계산 (옵션)
        weight_NT = None
        try:
            geo_path  = os.path.join(BASE_DIR, "geometry_temp.json")
            frame_path= os.path.join(BASE_DIR, "frame.json")
            if os.path.isfile(geo_path) and os.path.isfile(frame_path):
                with open(geo_path, "r", encoding="utf-8") as f:  geo = json.load(f)
                with open(frame_path, "r", encoding="utf-8") as f: frm = json.load(f)
                w_circ = float(geo.get("w_circumference", 0.0))
                wpm    = float(frm.get("weight_per_metter", 0.0))
                NT_frame = sum(1 for g in edge_group_name if g == "bottom ring")
                if wpm > 0.0 and NT_frame > 0 and w_circ > 0.0:
                    weight_NT = wpm * w_circ / NT_frame
        except Exception:
            pass

        meshinfo = {
            "Lines_pipe_top": lines_pipe_top,
            "numberOfLines_pipe_top": n_pipe,
            "Lines_netting": lines_netting,
            "numberOfLines_netting": n_twine,
            "surfs_netting": surfs_from_net,              # [[n1,n2,n3], [n1,n2,n3,n4], ...]
            "numberOfsurfs_netting": n_surfs,
            "Nodes": nodes_xyz
        }
        if weight_NT is not None:
            meshinfo["weight_NT"] = weight_NT  # Code_Aster: sinkF1에서 FZ=-meshinfo['weight_NT']

        os.makedirs(os.path.join(BASE_DIR, "saves"), exist_ok=True)
        with open(os.path.join(BASE_DIR, "saves", "meshinfo.json"), "w", encoding="utf-8") as f:
            json.dump(meshinfo, f, indent=2, ensure_ascii=False)
        print("✅ meshinfo.json saved.")
    except Exception as e:
        print(f"⚠️ meshinfo.json skipped: {e}")

    print(f"✅ VTU exported: {os.path.abspath(outpath)}")
    return outpath

if __name__ == "__main__":
    export_all_to_vtu()
