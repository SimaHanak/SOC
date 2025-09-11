import matplotlib.pyplot as plt
import csv
import numpy as np

with open(r"C:/Users/shana/Documents/[01] Studium/SOČ/trajectory.csv", 'r') as file:
#with open(r"/home/shanak/Documents/[01] Studium/SOČ/trajectory.csv", 'r') as file:
    reader = csv.reader(file)
    data = list(reader)
    data = [list(map(float, row)) for row in data[1:-1] if len(row) == 8]
    data = np.array(data)

print("Plotting", len(data), "points")

fig = plt.figure(figsize=(10, 6))
ax1 = fig.add_subplot(121, projection='3d')
ax1.set_xlabel('r*sin(phi)')
ax1.set_ylabel('r*cos(phi)')
ax1.set_zlabel('z')
ax1.set_title('3D Trajectory')

ax2 = fig.add_subplot(122)
ax2.set_xlabel('r')
ax2.set_ylabel('ur')
ax2.set_title('Poincare map')

# ax1.scatter(data[:, 2]*np.sin(data[:, 1]), data[:, 2]*np.cos(data[:, 1]), data[:, 3], s=1, c=data[:, 0], cmap='viridis', marker='o')
# poincare_map = []
# margin = 1e-2
# for i in range(len(data)):
#     if - margin < data[i, 3]%np.pi < margin:
#         poincare_map.append(data[i])
# poincare_map = np.array(poincare_map)
poincare_map = data
ax2.scatter(poincare_map[:, 2], poincare_map[:, 6], s=1, c=poincare_map[:, 0], cmap='viridis', marker='o')
# plt.colorbar(ax1.collections[0], label='t')
plt.show()
