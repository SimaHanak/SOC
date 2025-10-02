import matplotlib.pyplot as plt
import csv
import numpy as np
import os
from matplotlib.animation import FuncAnimation

list_dir = os.listdir("C:\\Users\\simon\\Documents\\[01] School\\[02] SOC\\SOC")
for i, dir in enumerate(list_dir):
    print(f"{i}: {dir}")

folder_number = int(input("Which folder you want to use? "))

with open(f"C:\\Users\\simon\\Documents\\[01] School\\[02] SOC\\SOC\\{list_dir[folder_number]}\\trajectory.csv", 'r') as tra_file:
#with open(r"/home/shanak/Documents/[01] Studium/SOČ/trajectory.csv", 'r') as file:
    reader = csv.reader(tra_file)
    data = list(reader)
    data = [list(map(float, row)) for row in data[1:-1] if len(row) == 8]
    data = np.array(data)

with open(f"C:\\Users\\simon\\Documents\\[01] School\\[02] SOC\\SOC\\{list_dir[folder_number]}\\rotation_numbers.csv", 'r') as rot_file:
    reader = csv.reader(rot_file)
    rotation_numbers_raw = list(reader)
    rotation_numbers = [[]]
    index = 0
    for num in rotation_numbers_raw:
        if len(num) > 1:
            if rotation_numbers[index] != []:
                index += 1
                rotation_numbers.append([])
        elif len(num) == 0:
            continue
        else:
            rotation_numbers[index].append(float(num[0]))

if rotation_numbers[-1] == []:
    rotation_numbers = rotation_numbers[:-1]

rotation_numbers = np.array(rotation_numbers, dtype=float)
for i in range(len(rotation_numbers)):
    if len(rotation_numbers[i]) > 0:
        rotation_numbers[i] = np.average(rotation_numbers[i])

print("Plotting", len(data), "points")

# fig = plt.figure(figsize=(10, 6))
# ax1 = fig.add_subplot(121)
# ax1.set_xlabel('x')
# ax1.set_ylabel('y')


# ax2 = fig.add_subplot(122)
# ax2.set_xlabel('r')
# ax2.set_ylabel('ur')
# ax2.set_title('Poincare map')

# # ax1.scatter(data[:, 2]*np.sin(data[:, 1]), data[:, 2]*np.cos(data[:, 1]), data[:, 3], s=1, c=data[:, 0], cmap='viridis', marker='o')
# # poincare_map = []
# # margin = 1e-2
# # for i in range(len(data)):
# #     if - margin < data[i, 3]%np.pi < margin:
# #         poincare_map.append(data[i])
# # poincare_map = np.array(poincare_map)
# ax1.plot(rotation_numbers, 'o-')

poincare_map = data
# ax2.scatter(poincare_map[:, 2], poincare_map[:, 6], s=0.1, cmap='viridis', marker='o') # c=poincare_map[:, 0]
# # plt.colorbar(ax1.collections[0], label='t')
# plt.show()


fig, ax = plt.subplots()
scat = ax.scatter([], [], s=1, cmap='viridis', marker='o')

ax.set_xlim(np.min(poincare_map[:, 2]), np.max(poincare_map[:, 2]))
ax.set_ylim(np.min(poincare_map[:, 6]), np.max(poincare_map[:, 6]))
ax.set_xlabel("r (col 2)")
ax.set_ylabel("ur (col 6)")

# update function for animation
def update(frame):
    x = poincare_map[:frame, 2]
    y = poincare_map[:frame, 6]
    scat.set_offsets(np.column_stack((x, y)))
    return scat,

# create animation
ani = FuncAnimation(fig, update, frames=len(poincare_map), interval=5, blit=True)

plt.show()
