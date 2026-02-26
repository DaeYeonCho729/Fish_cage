import numpy as np
import matplotlib.pyplot as plt

point1 = np.array((0, 2, 0))
point2 = np.array((2, 0, 0))
fluid_velocity = np.array((1,0, 0))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot points
ax.scatter(point1[0], point1[1], point1[2])
ax.scatter(point2[0], point2[1], point2[2])

# Element vector
elem_vec = point2 - point1
ax.quiver(point1[0], point1[1], point1[2],
          elem_vec[0], elem_vec[1], elem_vec[2])

# Fluid velocity vector
ax.quiver(0, 0, 0,
          fluid_velocity[0], fluid_velocity[1], fluid_velocity[2])

# Labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_zlim(-10, 10)

plt.title("3D visualization of line element and fluid velocity")
plt.show()