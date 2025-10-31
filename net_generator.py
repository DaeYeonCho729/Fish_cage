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

    with open("Fish_Cage/bottom_net.json", "w", encoding="utf-8") as f:
        json.dump({"points": bottom_net_pts}, f, indent=2, ensure_ascii=False)

    # 전체 점/엣지 저장: net_points.json
    with open("Fish_Cage/net_points.json", "w", encoding="utf-8") as f:
        json.dump({
            "points": points, "edges": edges,
            "z_float": float(z_float), "z_weight": float(z_weight)
        }, f, indent=2, ensure_ascii=False)

    print("✅ side_net.json / bottom_net.json (좌표만) / net_points.json saved")
    return points, edges, groups