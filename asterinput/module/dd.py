import numpy as np
from numpy import pi
point1 = np.array((0, 2, 0))
point2 = np.array((2, 0, 0))
fluid_velocity = np.array((1,0, 0))

a = point2 - point1
L = np.linalg.norm(a)
U = np.linalg.norm(fluid_velocity)

unit_normal_vector = a / (L)

if np.dot(unit_normal_vector, fluid_velocity) < 0:
    unit_normal_vector = - unit_normal_vector

coin_alpha = np.dot(unit_normal_vector, fluid_velocity) / (U)
inflow_angle = np.arccos(coin_alpha) # 도 -> 라디안
attack_angle = (pi/2) - inflow_angle
attack_angle_deg = np.degrees(attack_angle)

print("attack_angle (rad):", attack_angle)
print("attack_angle (deg):", attack_angle_deg)
print("inflow_angle (rad):", inflow_angle)
print("inflow_angle (deg):", np.degrees(inflow_angle))