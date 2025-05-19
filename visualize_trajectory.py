import matplotlib.pyplot as plt
import csv
import numpy as np

with open(r"C:\Users\simon\Documents\[01] School\[02] SOC\SOC\trajectory.csv", 'r') as file:
    reader = csv.reader(file)
    data = list(reader)
    data = [list(map(float, row)) for row in data]
    data = np.array(data)

print(data[0])

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('r*sin(phi)')
ax.set_ylabel('r*cos(phi)')
ax.set_zlabel('z')
ax.set_title('3D Trajectory')

ax.scatter(data[:, 2]*np.sin(data[:, 1]), data[:, 2]*np.cos(data[:, 1]), data[:, 3], c=data[:, 0], cmap='viridis', marker='o')
plt.colorbar(ax.collections[0], label='t')
plt.show()
