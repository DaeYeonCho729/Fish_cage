import numpy as np
import math

def generate_mooring_frame(sides, circumference, z_height, nodes_per_side, dedup_tol=None):
    """
    - 각 변마다 '코너 포함' 분할(모든 변에서 j=0 포함) → 코너-첫 내부점 edge가 누락되지 않음
    - 좌표 dedup은 나중에 수행(코너는 병합됨)
    - dedup_tol이 None이면, 최소 세그먼트 길이의 1e-8 배로 자동 설정(인접점이 합쳐지지 않도록 매우 작게)
    """
    # 안전장치
    sides = int(sides)
    nodes_per_side = max(2, int(nodes_per_side))  # 2 미만이면 edge가 안 생김

    R = circumference / (2.0 * math.pi)
    angle_step = 2.0 * math.pi / sides

    # 꼭짓점(코너)
    corner = [np.array([R*math.cos(i*angle_step), R*math.sin(i*angle_step), float(z_height)]) for i in range(sides)]

    raw_points = []
    edges = []
    start = 0

    # === 모든 변에서 j=0 포함(코너 포함) ===
    for i in range(sides):
        p0 = corner[i]
        p1 = corner[(i + 1) % sides]

        for j in range(0, nodes_per_side):
            t = j / (nodes_per_side - 1)
            p = (1.0 - t) * p0 + t * p1
            raw_points.append([float(p[0]), float(p[1]), float(p[2])])

        end = len(raw_points)
        # 이 변의 연속 edge
        for k in range(start, end - 1):
            edges.append((k, k + 1))
        start = end

    # ※ 별도의 '마지막→첫 점' 닫기 edge는 불필요(각 변이 코너까지 이미 연결)

    # ==== 좌표 중복 병합 ====
    # 자동 tol: 최소 세그먼트 길이 기반 (아주 작게)
    if dedup_tol is None:
        chord = 2.0 * R * math.sin(math.pi / sides)  # 한 변의 현 길이
        min_seg = chord / (nodes_per_side - 1)
        dedup_tol = max(1e-12, min_seg * 1e-8)

    inv = 1.0 / dedup_tol
    def qkey(p):
        return (round(p[0] * inv), round(p[1] * inv), round(p[2] * inv))

    key_to_new = {}
    new_points = []
    old_to_new = [None] * len(raw_points)

    for i, p in enumerate(raw_points):
        key = qkey(p)
        if key in key_to_new:
            old_to_new[i] = key_to_new[key]
        else:
            new_idx = len(new_points)
            key_to_new[key] = new_idx
            new_points.append(p)
            old_to_new[i] = new_idx

    # 엣지 리맵 + 자기자신/중복 제거
    seen = set()
    unique_edges = []
    for a, b in edges:
        a2 = old_to_new[a]; b2 = old_to_new[b]
        if a2 == b2:
            continue
        key = (a2, b2)
        if key in seen:
            continue
        seen.add(key)
        unique_edges.append(key)

    print(f"✅ Mooring Frame: {len(new_points)} unique nodes (from {len(raw_points)} raw), {len(unique_edges)} edges")
    return new_points, unique_edges
