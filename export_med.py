import os, json, glob
import salome
salome.salome_init()

import SMESH
from salome.smesh import smeshBuilder

BASE_DIR  = "Fish_Cage"
LINKS_DIR = os.path.join(BASE_DIR, "links")

def jload(path, default=None):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {} if default is None else default

def export_all_to_med(outpath="fish_cage.med"):
    # --- SALOME mesh context ---
    smesh = smeshBuilder.New()
    mesh  = smesh.Mesh()

    # --- global registries ---
    node_ids:   list[int]   = []     # Salome node IDs (global)
    node_coords: list[list] = []     # 각 글로벌 노드의 [x,y,z]
    edge_ids:   list[int]   = []     # Salome edge IDs (global)
    edge_pairs: list[tuple] = []     # global edge index -> (min(gi,gj), max(gi,gj))

    node_groups: dict[str, set|list] = {}  # {label: set(global_node_index)}
    edge_groups: dict[str, list]     = {}  # {label: list(global_edge_index)}

    # --- helpers ---
    def _map_idx(idx, section_pts_len, baseN):
        """
        섹션 로컬(0..section_pts_len-1)이면 baseN+idx,
        그 외(음수/범위 밖)는 '이미 존재하는 글로벌 인덱스'로 간주하여 그대로 사용.
        """
        i = int(idx)
        return (baseN + i) if (0 <= i < section_pts_len) else i

    def _add_section(points, edges, node_group_map=None, edge_group_map=None):
        """
        points: list[[x,y,z]]           -> 새 글로벌 노드로 추가
        edges:  list[[i,j]] (혼합 인덱스 허용: 섹션 로컬/글로벌)
        node_group_map: {label: [로컬/글로벌 혼합 노드 인덱스]}
        edge_group_map: {label: [섹션 로컬 엣지 인덱스(0..)]}
        """
        baseN = len(node_ids)
        secN  = len(points or [])
        baseE = len(edge_ids)

        # 1) add points (+좌표 저장)
        for p in (points or []):
            x, y, z = float(p[0]), float(p[1]), float(p[2])
            nid = mesh.AddNode(x, y, z)
            node_ids.append(nid)
            node_coords.append([x, y, z])

        # 2) add edges (혼합 인덱스 지원) + edge_pairs 기록
        for e in (edges or []):
            i, j = int(e[0]), int(e[1])
            gi = _map_idx(i, secN, baseN)
            gj = _map_idx(j, secN, baseN)
            if gi < 0 or gj < 0 or gi >= len(node_ids) or gj >= len(node_ids):
                continue
            eid = mesh.AddEdge([node_ids[gi], node_ids[gj]])
            edge_ids.append(eid)
            edge_pairs.append(tuple(sorted((gi, gj))))

        # 3) node groups (혼합 인덱스 지원)
        if node_group_map:
            for label, locs in node_group_map.items():
                if not locs:
                    continue
                G = node_groups.setdefault(label, set())
                for li in locs:
                    gi = _map_idx(int(li), secN, baseN)
                    if 0 <= gi < len(node_ids):
                        G.add(gi)

        # 4) edge groups (이 섹션에서 추가된 엣지의 로컬 인덱스 → 글로벌 인덱스)
        if edge_group_map:
            for label, locs in edge_group_map.items():
                if not locs:
                    continue
                L = edge_groups.setdefault(label, [])
                for li in locs:
                    ge = baseE + int(li)
                    if 0 <= ge < len(edge_ids):
                        L.append(ge)

    # ---------------- Sections ----------------
    # 1) Net  (z==0 → floater node, else net node; floater↔net 엣지는 net edge)
    net = jload(os.path.join(BASE_DIR, "net_points.json"), {})
    net_pts   = net.get("points", []) or []
    net_edges = net.get("edges", [])  or []

    baseN = len(node_ids)

    # 1-1) 노드 추가 + z기준 분류 (글로벌 인덱스 집합으로 저장)
    floater_nodes_from_net, net_nodes_from_net = set(), set()
    for i, p in enumerate(net_pts):
        x, y, z = float(p[0]), float(p[1]), float(p[2])
        nid = mesh.AddNode(x, y, z)
        node_ids.append(nid)
        node_coords.append([x, y, z])
        gi = baseN + i
        if abs(z) < 1e-8:
            floater_nodes_from_net.add(gi)   # z==0 → floater node
        else:
            net_nodes_from_net.add(gi)       # z!=0 → net node

    # 1-2) 엣지 추가(로컬/글로벌 혼재 안전 매핑) + 분류
    def _map_net_idx(idx: int) -> int:
        i = int(idx)
        return (baseN + i) if (0 <= i < len(net_pts)) else i

    floater_edges_from_net, net_edges_from_net = [], []
    for e in net_edges:
        i, j    = int(e[0]), int(e[1])
        gi, gj  = _map_net_idx(i), _map_net_idx(j)
        if gi < 0 or gj < 0 or gi >= len(node_ids) or gj >= len(node_ids):
            continue
        eid = mesh.AddEdge([node_ids[gi], node_ids[gj]])
        edge_ids.append(eid)
        ge = len(edge_ids) - 1
        edge_pairs.append(tuple(sorted((gi, gj))))

        # 규칙: 두 끝이 모두 floater-node일 때만 floater 엣지, 나머지는 전부 net 엣지
        if (gi in floater_nodes_from_net) and (gj in floater_nodes_from_net):
            floater_edges_from_net.append(ge)
        else:
            net_edges_from_net.append(ge)

    # 1-3) 그룹 반영 (기존 그룹과 합집합)
    if floater_nodes_from_net:
        node_groups["floater node"] = node_groups.get("floater node", set()) | floater_nodes_from_net
    if net_nodes_from_net:
        node_groups["net node"]     = node_groups.get("net node", set())     | net_nodes_from_net
    if floater_edges_from_net:
        edge_groups["floater"]      = edge_groups.get("floater", [])         + floater_edges_from_net
    if net_edges_from_net:
        edge_groups["net twine"]    = edge_groups.get("net twine", [])       + net_edges_from_net

    # 2) Floating collar (optional)
    flo = jload(os.path.join(BASE_DIR, "float_temp.json"), {})
    flo_pts   = flo.get("points", [])
    flo_edges = flo.get("edges", [])
    if flo_pts or flo_edges:
        node_map = {"floater node": list(range(len(flo_pts)))} if flo_pts else None
        edge_map = {"floater": list(range(len(flo_edges)))} if flo_edges else None
        _add_section(flo_pts, flo_edges, node_map, edge_map)

    # 3) Bottom ring / collar (optional) — 저장 시점 baseN을 현재 장면으로 시프트
    bot = jload(os.path.join(BASE_DIR, "bottom_collar_temp.json"), {})
    bot_pts   = bot.get("points_new", []) or bot.get("points", [])
    bot_edges = bot.get("ring_edges", []) or bot.get("edges", [])
    bot_ring_indices = bot.get("ring_indices", [])
    bot_baseN = bot.get("baseN", None)

    if bot_pts or bot_edges or bot_ring_indices:
        curr_before = len(node_ids)
        def _shift(n: int) -> int:
            if bot_baseN is None:
                return int(n)
            n = int(n)
            shift = curr_before - int(bot_baseN)
            return n + shift if n >= int(bot_baseN) else n

        bot_edges_shifted = [[_shift(e[0]), _shift(e[1])] for e in (bot_edges or [])]
        bot_ring_shifted  = [_shift(i) for i in (bot_ring_indices or [])]

        node_map = {"bottom ring node": bot_ring_shifted} if bot_ring_shifted else None
        edge_map = {"bottom ring": list(range(len(bot_edges_shifted)))} if bot_edges_shifted else None
        _add_section(bot_pts, bot_edges_shifted, node_map, edge_map)

    # 4) Mooring frame (optional) — 저장 시점 baseN을 현재 장면으로 시프트
    moor = jload(os.path.join(BASE_DIR, "mooring_temp.json"), {})
    moor_pts   = moor.get("points", [])
    moor_edges = moor.get("edges", [])
    moor_baseN = moor.get("baseN", None)

    if moor_pts or moor_edges:
        curr_before = len(node_ids)
        def _shift_moor(n: int) -> int:
            n = int(n)
            if moor_baseN is None:
                return n
            shift = curr_before - int(moor_baseN)
            return n + shift if n >= int(moor_baseN) else n

        moor_edges_shifted = [[_shift_moor(e[0]), _shift_moor(e[1])] for e in (moor_edges or [])]
        node_map = {"mooring frame node": list(range(len(moor_pts)))} if moor_pts else None
        edge_map = {"mooring frame": list(range(len(moor_edges_shifted)))} if moor_edges_shifted else None
        _add_section(moor_pts, moor_edges_shifted, node_map, edge_map)

    # 5) Links (optional) — 저장 시점 baseN을 현재 장면으로 시프트 + buoy/anchor 노드 그룹 채움
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

    # 필요 없는 노드그룹 제거: side rope node, anchor line node는 만들지 않음
    link_nodes = {
        "buoy1 node": set(),
        "buoy2 node": set(),
        "mooring line node": set(),
        # anchor node는 최종에 폴리라인별로 1개씩 계산
    }

    # 앵커 라인 폴리라인 후보 저장용 (A/B를 분리해 페어링)
    _anchor_lines_A = []  # list[ (ns_shifted, pts_len, baseN_before) ]
    _anchor_lines_B = []  # list[ (ns_shifted, pts_len, baseN_before) ]

    if os.path.isdir(LINKS_DIR):
        for f in sorted(glob.glob(os.path.join(LINKS_DIR, "*.json"))):
            data = jload(f, {})
            pts  = data.get("extra_points", []) or []
            polylines = data.get("polylines", []) or []
            saved_baseN = data.get("baseN", None)

            curr_before = len(node_ids)   # 이 파일 시작 시점의 글로벌 노드 수

            def _S(n: int) -> int:
                n = int(n)
                if saved_baseN is None:
                    return n  # 이미 글로벌 인덱스로 저장된 경우
                shift = curr_before - int(saved_baseN)
                return n + shift if n >= int(saved_baseN) else n

            edges_local = []
            edge_group_local = {}

            for pl in polylines:
                ns_raw = [int(x) for x in (pl.get("nodes", []) or [])]
                if len(ns_raw) < 2:
                    continue

                ns = [_S(x) for x in ns_raw]   # 모든 노드에 시프트 적용

                # 엣지 생성(연속)
                start_e = len(edges_local)
                for i in range(len(ns) - 1):
                    edges_local.append([ns[i], ns[i + 1]])

                # 그룹명 매핑
                label = key2edge.get(pl.get("group", ""), None)
                if label:
                    edge_idx_list = list(range(start_e, len(edges_local)))
                    edge_group_local.setdefault(label, []).extend(edge_idx_list)

                # ---- 링크 노드 그룹 수집(필요한 것만) ----
                g = pl.get("group", "")
                if g == "BUOY_TETHER":
                    link_nodes["buoy1 node"].add(ns[-1])
                elif g == "BUOY_LINE":
                    link_nodes["buoy2 node"].add(ns[-1])
                elif g == "DISTANCE_ROPE":
                    link_nodes["mooring line node"].update(ns)

                # ★ 앵커 라인: A/B를 분리해서 폴리라인 단위로 저장
                if g == "ANCHOR_LINE_A":
                    _anchor_lines_A.append((ns, len(pts), curr_before))
                elif g == "ANCHOR_LINE_B":
                    _anchor_lines_B.append((ns, len(pts), curr_before))

            # extra_points는 로컬 0..N-1 새 노드로 간주 → 여기서 실제 추가/매핑 수행
            _add_section(pts, edges_local, None, edge_group_local if edge_group_local else None)

    # ---- anchor node 최종 결정: 각 A/B "페어"마다 z 최저 1개씩 선택 ----
    def _to_global(ns, pts_len, baseN_before):
        out = []
        for i in ns:
            ii = int(i)
            if 0 <= ii < pts_len:
                out.append(baseN_before + ii)
            else:
                out.append(ii)
        return out

    def _min_z_node(indices):
        best, best_z = None, None
        for idx in indices:
            i = int(idx)
            if 0 <= i < len(node_coords):
                z = node_coords[i][2]
                if (best_z is None) or (z < best_z):
                    best_z, best = z, i
        return best

    # A/B 각각 폴리라인들을 글로벌 인덱스로 변환
    A_gl = [_to_global(ns, pts_len, baseN_before) for (ns, pts_len, baseN_before) in _anchor_lines_A]
    B_gl = [_to_global(ns, pts_len, baseN_before) for (ns, pts_len, baseN_before) in _anchor_lines_B]

    anchor_final = set()
    max_pairs = max(len(A_gl), len(B_gl))
    for k in range(max_pairs):
        cand = set()
        if k < len(A_gl): cand.update(A_gl[k])
        if k < len(B_gl): cand.update(B_gl[k])
        winner = _min_z_node(cand)
        if winner is not None:
            anchor_final.add(winner)

    if anchor_final:
        node_groups["anchor node"] = anchor_final
    else:
        if "anchor node" in node_groups:
            del node_groups["anchor node"]

    # ---- 링크 기반 나머지 노드 그룹 반영 ----
    for label, S in link_nodes.items():
        if not S:
            continue
        G = node_groups.setdefault(label, set())
        G.update(i for i in S if 0 <= int(i) < len(node_ids))

    # ===================== DEDUP: bottom ring vs net twine =====================
    def _dedup_bottom_vs_net(edge_pairs, edge_groups, node_groups,
                             net_name="net twine", bottom_name="bottom ring"):
        """
        edge_pairs: list[tuple(int,int)]   # 각 글로벌 엣지의 (gi, gj)
        edge_groups: dict[str, list[int]]  # 그룹명 -> 글로벌 엣지 인덱스
        node_groups: dict[str, list|set]   # 그룹명 -> 글로벌 노드 인덱스
        """
        if bottom_name not in edge_groups or net_name not in edge_groups:
            return

        # 1) bottom ring이 보유한 엣지쌍 집합
        bottom_set = set()
        for eid in edge_groups.get(bottom_name, []):
            if 0 <= eid < len(edge_pairs):
                bottom_set.add(edge_pairs[eid])

        # 2) net twine 엣지에서 bottom과 겹치는 엣지 제거
        old_net_eids = [eid for eid in edge_groups.get(net_name, []) if 0 <= eid < len(edge_pairs)]
        new_net_eids = [eid for eid in old_net_eids if edge_pairs[eid] not in bottom_set]
        edge_groups[net_name] = new_net_eids

        # 3) 남은 net twine 엣지로부터 노드 재계산 → net twine 노드 그룹 갱신
        keep_nodes = set()
        for eid in new_net_eids:
            gi, gj = edge_pairs[eid]
            keep_nodes.add(gi); keep_nodes.add(gj)

        if net_name in node_groups:
            orig = node_groups[net_name]
            if isinstance(orig, set):
                node_groups[net_name] = set([nid for nid in orig if nid in keep_nodes])
            else:
                node_groups[net_name] = sorted([nid for nid in orig if nid in keep_nodes])

        # 4) bottom ring 노드 그룹 보정: bottom 엣지에서 직접 수집
        if bottom_name in node_groups:
            bottom_nodes = set()
            for eid in edge_groups.get(bottom_name, []):
                if 0 <= eid < len(edge_pairs):
                    gi, gj = edge_pairs[eid]
                    bottom_nodes.add(gi); bottom_nodes.add(gj)
            if isinstance(node_groups[bottom_name], set):
                node_groups[bottom_name] = bottom_nodes
            else:
                node_groups[bottom_name] = sorted(bottom_nodes)

    _dedup_bottom_vs_net(edge_pairs, edge_groups, node_groups)

    # all node
    node_groups["all node"] = set(range(len(node_ids)))

    # --- SMESH 그룹 생성 ---
    for label, idxs in node_groups.items():
        idxs_sorted = sorted(list(idxs)) if isinstance(idxs, set) else sorted(idxs)
        if not idxs_sorted:
            continue
        g = mesh.CreateEmptyGroup(SMESH.NODE, label)
        g.Add([node_ids[i] for i in idxs_sorted])

    for label, eidxs in edge_groups.items():
        if not eidxs:
            continue
        g = mesh.CreateEmptyGroup(SMESH.EDGE, label)
        g.Add([edge_ids[i] for i in eidxs])

    # --- MED 내보내기 ---
    outpath = os.path.abspath(outpath)
    mesh.Compute()
    mesh.ExportMED(outpath, 0)
    print(f"✅ MED exported: {outpath}")

        # ---------------- meshinfo.json (MED 노드 번호로 저장) ----------------
    try:
        # 0) 준비
        os.makedirs(os.path.join(BASE_DIR, "saves"), exist_ok=True)

        # net_points.json에서 z_float, surfs_netting(0-베이스) 읽기
        net_json_path = os.path.join(BASE_DIR, "net_points.json")
        net_raw = jload(net_json_path, {})
        z_float = float(net_raw.get("z_float", 0.0))
        surfs_from_net = net_raw.get("surfs_netting", []) or []

        # 1) MED 노드 번호(ID) ↔ 좌표
        #    - node_ids: 우리가 AddNode할 때 받은 SMESH Node ID(보통 1..N)
        #    - coords_by_med[nid] = (x,y,z)
        coords_by_med = {}
        for nid in node_ids:
            x, y, z = mesh.GetNodeXYZ(nid)
            coords_by_med[int(nid)] = [float(x), float(y), float(z)]

        # 2) edge 그룹을 MED 노드 번호로 변환
        #    edge_ids[global_edge_index] -> eid(SMESH Elem ID for EDGE)
        #    mesh.GetElemNodes(eid) -> [nid1, nid2]
        def _edge_list_for_group(gname: str):
            e_globs = edge_groups.get(gname, []) or []
            out = []
            for ge in e_globs:
                if 0 <= ge < len(edge_ids):
                    eid = edge_ids[ge]
                    nids = mesh.GetElemNodes(eid)  # [nid1, nid2] (1-베이스)
                    if len(nids) == 2:
                        out.append([int(nids[0]), int(nids[1])])
            return out

        lines_floater_all = _edge_list_for_group("floater")
        lines_netting_all = _edge_list_for_group("net twine")
        lines_braket_all  = _edge_list_for_group("bracket")   # ← bracket 그룹을 Line_braket로

        # 3) Lines_pipe_top: floater 중 z≈z_float인 것만
        def _approx(a, b, tol=1e-6): return abs(a - b) <= tol
        lines_pipe_top = []
        for a_med, b_med in lines_floater_all:
            za = coords_by_med[a_med][2]; zb = coords_by_med[b_med][2]
            if _approx(za, z_float) and _approx(zb, z_float):
                lines_pipe_top.append([a_med, b_med])

        # 4) surfs_netting(0-베이스 글로벌 인덱스) → MED 노드 번호로 매핑
        #    주의: net_points.json의 points 추가가 export_med의 노드 추가와 동일 순서이므로
        #          "글로벌 0-베이스 인덱스 i" -> "node_ids[i] (MED Node ID)" 로 매핑 가능.
        surfs_med = []
        for face in surfs_from_net:
            # face는 [i,j,k,(l)] 형태의 0-베이스 글로벌 인덱스
            try:
                surfs_med.append([int(node_ids[int(i)]) for i in face])
            except Exception:
                # 범위를 벗어나는 인덱스가 있으면 스킵
                pass

        # 5) Nodes: MED 노드 번호 순서로 좌표 정리(키 정렬)
        med_ids_sorted = sorted(coords_by_med.keys())
        nodes_xyz_by_med = [[coords_by_med[n][0], coords_by_med[n][1], coords_by_med[n][2]]
                            for n in med_ids_sorted]

        # 6) meshinfo 구성 및 저장
        meshinfo = {
            "Lines_pipe_top": lines_pipe_top,
            "numberOfLines_pipe_top": len(lines_pipe_top),

            "Lines_netting": lines_netting_all,
            "numberOfLines_netting": len(lines_netting_all),

            "Line_braket": lines_braket_all,                    
            "numberOfLine_braket": len(lines_braket_all),       

            "surfs_netting": surfs_med,
            "numberOfsurfs_netting": len(surfs_med),

            "Nodes": nodes_xyz_by_med,
            "MED_node_ids": med_ids_sorted  # ← 원한다면 남겨두면 디버깅에 유용
        }

        out_meshinfo = os.path.join(BASE_DIR, "saves", "meshinfo.json")
        with open(out_meshinfo, "w", encoding="utf-8") as f:
            json.dump(meshinfo, f, indent=2, ensure_ascii=False)

        print(f"✅ meshinfo.json (MED numbering) saved: {out_meshinfo}")
    except Exception as e:
        print(f"⚠️ meshinfo.json (MED numbering) skipped: {e}")

    return outpath

if __name__ == "__main__":
    export_all_to_med("fish_cage.med")
