import numpy as np
import json, os
import vtk

def generate_net_nodes_and_edges(
    sides, f_circumference, w_circumference,
    z_float, z_weight, half_mesh_size,
    bottomnet_angle_deg  # ✅ 추가
):
    H = abs(z_float - z_weight)
    NN = int(np.ceil(H / half_mesh_size))
    z_levels = np.linspace(z_float, z_weight, NN + 1)

    NT_top = int(np.ceil(f_circumference / half_mesh_size))
    NT_bottom = int(np.ceil(w_circumference / half_mesh_size))
    NT_max = max(NT_top, NT_bottom, sides)

    segments_per_side = NT_max // sides
    if NT_max % sides != 0:
        segments_per_side += 1
        NT_max = segments_per_side * sides

    points = []
    edges = []
    prev_ring = []

    groups = {"top_ring": [], "bottom_ring": [], "bottom_plate": []}
    tol = 1e-6
    first_ring_saved = None
    
    for z in z_levels:
        C = f_circumference + (w_circumference - f_circumference) * ((z - z_float) / (z_weight - z_float))
        R = C / (2 * np.pi)

        vertices = []
        for i in range(sides):
            angle = 2 * np.pi * i / sides
            x = R * np.cos(angle)
            y = R * np.sin(angle)
            vertices.append([x, y])

        ring = []
        for i in range(sides):
            start = np.array(vertices[i])
            end = np.array(vertices[(i + 1) % sides])
            for j in range(segments_per_side):
                t = j / segments_per_side
                p = (1 - t) * start + t * end
                points.append([p[0], p[1], z])
                ring.append(len(points) - 1)

        NT = len(ring)
        for i in range(NT):
            edges.append((ring[i], ring[(i + 1) % NT]))

        if prev_ring:
            for i in range(NT):
                edges.append((prev_ring[i], ring[i]))
        
        # 첫 링 저장 (top)
        if first_ring_saved is None:
            first_ring_saved = list(ring)

        prev_ring = ring

    # === 바닥 그물판 (중앙으로 기울기 반영) ===
    if prev_ring:
        bottom_ring = prev_ring
        NT = len(bottom_ring)

        center_xy = np.array([0.0, 0.0])
        bottom_coords = np.array([points[idx][:2] for idx in bottom_ring])

        lengths = [np.linalg.norm(p - center_xy) for p in bottom_coords]
        max_radius = max(lengths)

        tilt_rad = np.radians(bottomnet_angle_deg)
        z_center = z_weight - max_radius * np.tan(tilt_rad)

        center_point = [0.0, 0.0, z_center]
        points.append(center_point)
        center_idx = len(points) - 1

        radial_div = max(1, int(np.ceil(max_radius / half_mesh_size)))
        ring_nodes = []

        for i in range(NT):
            p_outer = bottom_coords[i]
            nodes = [center_idx]
            for r in range(1, radial_div + 1):
                ratio = r / radial_div
                p_xy = center_xy + (p_outer - center_xy) * ratio
                z = z_center + (z_weight - z_center) * ratio
                points.append([p_xy[0], p_xy[1], z])
                nodes.append(len(points) - 1)
            ring_nodes.append(nodes)

        for ray in ring_nodes:
            for i in range(len(ray) - 1):
                edges.append((ray[i], ray[i + 1]))

        for step in range(1, radial_div + 1):
            for i in range(NT):
                ray_curr = ring_nodes[i]
                ray_next = ring_nodes[(i + 1) % NT]
                edges.append((ray_curr[step], ray_next[step]))

    for step in range(1, radial_div + 1):
            for i in range(NT):
                ray_curr = ring_nodes[i]
                ray_next = ring_nodes[(i + 1) % NT]
                edges.append((ray_curr[step], ray_next[step]))

    dedup_points = []
    old_to_dedup = {}
    seen = {}

    def key3(p, nd=6):
        return (round(float(p[0]), nd), round(float(p[1]), nd), round(float(p[2]), nd))

    for old_idx, p in enumerate(points):
        k = key3(p)
        if k in seen:
            old_to_dedup[old_idx] = seen[k]
        else:
            dedup_idx = len(dedup_points)
            dedup_points.append(p)
            seen[k] = dedup_idx
            old_to_dedup[old_idx] = dedup_idx

    import math

    def sort_key(i):
        x, y, z = dedup_points[i]
        ang = math.atan2(y, x) % (2 * math.pi)  # 0~2π 각도
        return (-z, ang)  # z는 큰 값 먼저(-z), 각도는 CCW(작은→큰)

    order = sorted(range(len(dedup_points)), key=sort_key)

    dedup_to_final = {dedup_idx: final_idx for final_idx, dedup_idx in enumerate(order)}
    old_to_final = {old_idx: dedup_to_final[old_to_dedup[old_idx]] for old_idx in range(len(points))}

    points = [dedup_points[di] for di in order]

    def remap_edge(e):
        a, b = int(e[0]), int(e[1])
        return (old_to_final[a], old_to_final[b])

    edges = [remap_edge(e) for e in edges]


    print(f"✅ 총 노드 수: {len(points)}, 총 Edge 수: {len(edges)}")

    # === 그룹 지정 (return 직전) ===
    if first_ring_saved is not None:
        groups["top_ring"] = list(first_ring_saved)
    if prev_ring:
        groups["bottom_ring"] = list(prev_ring)
    # 바닥판 전체: z ≤ z_weight + tol
    groups["bottom_plate"] = [i for i, p in enumerate(points) if p[2] <= z_weight + tol]

    print(f"✅ 총 노드 수: {len(points)}, 총 Edge 수: {len(edges)}")

    os.makedirs("Fish_Cage", exist_ok=True)

    # side_net / bottom_net 분리 (좌표만 저장)
    n = len(points)
    top_ring = set(groups.get("top_ring", []))
    bottom_plate = set(groups.get("bottom_plate", []))
    side_net_pts = [points[i] for i in range(n) if i not in top_ring and i not in bottom_plate]
    bottom_net_pts = [points[i] for i in bottom_plate]

    with open("Fish_Cage/side_net.json", "w", encoding="utf-8") as f:
        json.dump({"points": side_net_pts}, f, indent=2, ensure_ascii=False)

    # === bottom_net.json: points + surfs(삼각/사각) 저장 ===
    # bottom_faces 는 최종(전역) 인덱스 기준 → bottom 영역 로컬 인덱스로 변환
    g2l = {}
    local_pts = []
    for gi in sorted(bottom_plate):
        g2l[gi] = len(local_pts)
        local_pts.append(points[gi])

    bottom_surfs_local = []
    # bottom_faces 가 없는 경우도 있으니 방어
    if 'bottom_faces' in locals() and bottom_faces:
        for face in bottom_faces:
            # face 가 전부 bottom_plate에 속하는 경우만 기록
            if all(v in bottom_plate for v in face):
                bottom_surfs_local.append([g2l[v] for v in face])

    with open("Fish_Cage/bottom_net.json", "w", encoding="utf-8") as f:
        json.dump(
            {"points": local_pts, "surfs": bottom_surfs_local},
            f, indent=2, ensure_ascii=False
        )

    # =========================
    # surfs_netting 생성 (층별 연속 인덱스 가정)
    #  - 링 수: NN+1
    #  - 각 링 노드 수: NT_max
    #  - 각 링의 노드는 [k*NT_max, (k+1)*NT_max) 범위로 저장되어 있다고 가정
    # =========================
    surfs = []

    # 방어적으로 필요한 값이 없으면 계산
    import math
    if 'NT_max' not in locals():
        NT_top    = int(math.ceil(f_circumference / half_mesh_size))
        NT_bottom = int(math.ceil(w_circumference / half_mesh_size))
        NT_max    = max(NT_top, NT_bottom, sides)

    if 'NN' not in locals():
        H  = abs(z_float - z_weight)
        NN = int(math.ceil(H / half_mesh_size))  # 수직 분할 수

    total_nodes = len(points)

    # 1) 링 인덱스 테이블 구성
    ring_indices = []
    for k in range(NN + 1):  # NN+1개의 링
        start = k * NT_max
        end   = start + NT_max
        if end > total_nodes:
            # 데이터가 덜 쌓였으면 중단
            break
        ring_indices.append(list(range(start, end)))

    # 2) 인접 링 간 사각 패널: (i, i+1, i+1, i)
    for k in range(len(ring_indices) - 1):
        R0 = ring_indices[k]
        R1 = ring_indices[k + 1]
        for i in range(NT_max):
            a = R0[i]
            b = R0[(i + 1) % NT_max]
            c = R1[(i + 1) % NT_max]
            d = R1[i]
            surfs.append([a, b, c, d])  # 사각

    # 3) 바닥 중앙 삼각(옵션): ring_nodes가 있으면 중심-첫 고리 삼각 추가
    #    ring_nodes = [[center, r1_0, r2_0, ...], ..., spoke별 리스트] 구조라면
    bottom_faces = []  # ← 나중에 bottom_net.json 저장에도 씀
    if 'ring_nodes' in locals() and ring_nodes:
        # old_to_final 없으면 항등 매핑
        if 'old_to_final' not in locals():
            old_to_final = {i: i for i in range(total_nodes)}

        # 방사선(스포크) → 최종 인덱스
        rays = [[old_to_final.get(idx, idx) for idx in ray] for ray in ring_nodes]

        # ----- (A) 삼각: 센터–첫 고리 -----
        center     = rays[0][0]
        first_ring = [ray[1] for ray in rays]
        NT_bot     = len(first_ring)
        for i in range(NT_bot):
            a = center
            b = first_ring[i]
            c = first_ring[(i + 1) % NT_bot]
            if len({a, b, c}) == 3:
                bottom_faces.append([a, b, c])

        # ----- (B) 사각: 고리 s ↔ s+1 -----
        # ring_nodes 구조상 rays[*][s] 가 안쪽, rays[*][s+1] 가 바깥
        # s=1..(radial_div-1) 범위가 "첫 고리~마지막 고리" 띠를 이룸
        for s in range(1, radial_div):
            for i in range(NT_bot):
                a = rays[i][s]
                b = rays[i][s + 1]
                c = rays[(i + 1) % NT_bot][s + 1]
                d = rays[(i + 1) % NT_bot][s]
                if len({a, b, c, d}) == 4:
                    bottom_faces.append([a, b, c, d])

    # 측면 사각(위에서 만든 surfs) + 바닥(삼각/사각) 모두 합치기
    surfs_netting = []
    for face in surfs:          # 측면 그물 패널마다
        surfs_netting.append([n + 1 for n in face] + ['s'])   # 노드 인덱스를 1증가 후 's' 태그 추가
    for face in bottom_faces:   # 바닥 그물 패널마다
        surfs_netting.append([n + 1 for n in face] + ['b'])   # 노드 인덱스를 1증가 후 'b' 태그 추가

    # net_points.json 파일에 저장 (points, edges 등과 함께)
    with open("Fish_Cage/net_points.json", "w", encoding="utf-8") as f:
        json.dump({
            "points": points,
            "edges": edges,
            "surfs_netting": surfs_netting,  # ★ 각 패널의 끝에 태그 포함
            "z_float": float(z_float),
            "z_weight": float(z_weight),
            "numerical_half_mesh_size": half_mesh_size  # ← 이 줄 추가
        }, f, indent=2, ensure_ascii=False)

    print("✅ side_net.json / bottom_net.json (좌표) / net_points.json(+surfs) saved")
    return points, edges, groups
