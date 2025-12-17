import numpy as np

def generate_cage_frame(sides, f_circumference, z_float, w_circumference, z_weight):
    f_radius = f_circumference / (2 * np.pi)
    w_radius = w_circumference / (2 * np.pi)
    angle_step = 2 * np.pi / sides

    float_coords = []
    weight_coords = []
    frame_edges = []

    for i in range(sides):
        angle = i * angle_step

        x_f = f_radius * np.cos(angle)
        y_f = f_radius * np.sin(angle)
        float_coords.append([x_f, y_f, z_float])

        x_w = w_radius * np.cos(angle)
        y_w = w_radius * np.sin(angle)
        weight_coords.append([x_w, y_w, z_weight])

        frame_edges.append((i, i + sides))  # 수직

    for i in range(sides):
        frame_edges.append((i, (i + 1) % sides))
        frame_edges.append((i + sides, ((i + 1) % sides) + sides))

    points = float_coords + weight_coords
    return points, frame_edges
