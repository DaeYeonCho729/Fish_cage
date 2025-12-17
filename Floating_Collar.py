import numpy as np
import json

def generate_outer_floater(sides, bracket_per_side, circumference, z_pos,
                           output_path="Fish_Cage/float_temp.json"):
    # 반지름 계산 (정다각형 둘레 = circumference)
    R = circumference / (2 * sides * np.sin(np.pi / sides))

    all_pts = []
    corner_idxs = []

    # 1) 꼭짓점 → 그 변의 bracket 점들 순으로 interleaved 생성
    for i in range(sides):
        θ0 = 2 * np.pi * i / sides
        θ1 = 2 * np.pi * (i + 1) / sides

        p0 = np.array([R * np.cos(θ0), R * np.sin(θ0)])
        p1 = np.array([R * np.cos(θ1), R * np.sin(θ1)])

        # 꼭짓점 추가
        corner_idxs.append(len(all_pts))
        all_pts.append([p0[0], p0[1], z_pos])

        # bracket_per_side 만큼 보간점 추가
        for j in range(bracket_per_side):
            t = (j + 1) / (bracket_per_side + 1)
            P = (1 - t) * p0 + t * p1
            all_pts.append([P[0], P[1], z_pos])

    total = len(all_pts)

    # 2) Skeleton edge: 꼭짓점끼리만 연결
    skeleton_edges = [
        (corner_idxs[i], corner_idxs[(i + 1) % sides])
        for i in range(sides)
    ]

    data = {
        "points": all_pts,
        "edges": skeleton_edges,
        "z_float": z_pos
    }
    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)

    print(f"✅ 외부 floater 저장: {sides}각형 + 변당 {bracket_per_side}개 보간 → 총 노드 {total}")
    return all_pts, skeleton_edges