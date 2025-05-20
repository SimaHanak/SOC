import matplotlib.pyplot as plt
import csv
import numpy as np

with open(r"/home/shanak/Documents/[01] Studium/SOČ/trajectory.csv", 'r') as file:
    reader = csv.reader(file)
    data = list(reader)
    data = [list(map(float, row)) for row in data]
    data = np.array(data)

print(data[0])

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

ax1.scatter(data[:, 2]*np.sin(data[:, 1]), data[:, 2]*np.cos(data[:, 1]), data[:, 3], c=data[:, 0], cmap='viridis', marker='o')
poincare_map = []
for i in range(len(data)):
    if round(data[i, 1], -1) == 0:
        poincare_map.append(data[i]) 
poincare_map = np.array(poincare_map)
ax2.scatter(poincare_map[:, 2], poincare_map[:, 6])
plt.colorbar(ax1.collections[0], label='t')
plt.show()
