import matplotlib.pyplot as plt
import csv
import numpy as np
import os
from matplotlib.animation import FuncAnimation

# list_dir = os.listdir("C:\\Users\\simon\\Documents\\[01] School\\[02] SOC\\SOC")
list_dir = os.listdir("/home/shanak/Documents/[01] Studium/SOČ/")
for i, dir in enumerate(list_dir):
    print(f"{i}: {dir}")

folder_number = int(input("Which folder you want to use? "))

if list_dir[folder_number].endswith(".csv"):
    file_path = f"/home/shanak/Documents/[01] Studium/SOČ/{list_dir[folder_number]}"
else:
    file_path = f"/home/shanak/Documents/[01] Studium/SOČ/{list_dir[folder_number]}/trajectory.csv"

#with open(f"C:\\Users\\simon\\Documents\\[01] School\\[02] SOC\\SOC\\{list_dir[folder_number]}\\trajectory.csv", 'r') as tra_file:
with open(file_path, 'r') as tra_file:
    reader = csv.reader(tra_file)
    data_raw = list(reader)
    data_raw.append([])
    data = []
    length = 0
    lengths = []
    header_r = []
    for row in data_raw[1:]:
        if len(row) == 2:
            data.append(list(map(float, row)))
            length += 1
        else:
            lengths.append(length)
            length = 0
    for row in data_raw:
        if len(row) == 9:
            header_r.append(float(row[2][2:]))
    data = np.array(data)


# with open(f"C:\\Users\\simon\\Documents\\[01] School\\[02] SOC\\SOC\\{list_dir[folder_number]}\\rotation_numbers.csv", 'r') as rot_file:
#     reader = csv.reader(rot_file)
#     rotation_numbers_raw = list(reader)
#     rotation_numbers_header = [float(head[2][2:]) for head in rotation_numbers_raw[::2]]
#     rotation_numbers_data = [float(num[0]) for num in rotation_numbers_raw[1::2]]
#     print({rotation_numbers_header[i]: rotation_numbers_data[i] for i in range(len(rotation_numbers_header))})

def compute_theta(r0, ur0, r, ur, rc):
    angle = np.arctan2(ur, r - rc) - np.arctan2(ur0, r0 - rc)
    return angle if angle >= 0 else angle + 2 * np.pi

def compute_center(r_values):
    r_min = np.min(r_values)
    r_max = np.max(r_values)
    return (r_min + r_max) / 2

center = compute_center(data[:lengths[0],0])
rotation_numbers = []
print("Center r =", center)
for j in range(len(lengths)):
    start_index = sum(lengths[:j])
    end_index = start_index + lengths[j]
    tmp_values = data[start_index:end_index]
    angle = 0
    for i in range(1, lengths[j]):
        angle += compute_theta(tmp_values[i-1,0], tmp_values[i-1,1], tmp_values[i,0], tmp_values[i,1], center)
    if lengths[j] != 0:
        rotation_numbers.append(angle / (2 * np.pi * (lengths[j]-1)))

print("Number of rotation numbers computed:", len(rotation_numbers))

print("Plotting", len(data), "points")

fig = plt.figure(figsize=(10, 6))
ax1 = fig.add_subplot(121)
ax1.set_xlabel('r')
ax1.set_ylabel('Rotation number')
ax1.tick_params(axis='x', labelrotation=45)

ax2 = fig.add_subplot(122)
ax2.set_xlabel('r')
ax2.set_ylabel('ur')
ax2.set_title('Poincare map')

# ax3 = fig.add_subplot(223)
# ax3.set_xlabel('r')
# ax3.set_ylabel('ur')
# ax3.set_title('Rotation numbers, python computed')  

if len(header_r) > len(rotation_numbers):
    header_r = header_r[:len(rotation_numbers)]
ax1.scatter(header_r, rotation_numbers, s=4, marker='o')

for i, rnum in enumerate(rotation_numbers):
    print(i, rnum)
ax2.scatter(data[:, 0], data[:, 1], s=0.1, cmap='viridis', marker='o') # c=poincare_map[:, 0]

# ax3.scatter(rotation_numbers_header, rotation_numbers, s=4, marker='o')
# plt.colorbar(ax1.collections[0], label='t')
plt.show()


# fig, ax = plt.subplots()
# scat = ax.scatter([], [], s=1, cmap='viridis', marker='o')

# ax.set_xlim(np.min(poincare_map[:, 2]), np.max(poincare_map[:, 2]))
# ax.set_ylim(np.min(poincare_map[:, 6]), np.max(poincare_map[:, 6]))
# ax.set_xlabel("r (col 2)")
# ax.set_ylabel("ur (col 6)")

# # update function for animation
# def update(frame):
#     x = poincare_map[:frame, 2]
#     y = poincare_map[:frame, 6]
#     scat.set_offsets(np.column_stack((x, y)))
#     return scat,

# # create animation
# ani = FuncAnimation(fig, update, frames=len(poincare_map), interval=5, blit=True)

# plt.show()
