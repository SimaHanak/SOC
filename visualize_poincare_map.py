import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from bokeh.plotting import figure, show, ColumnDataSource
from bokeh.layouts import row
from bokeh.models import ZoomInTool, ZoomOutTool

def load_data():
    list_dir = os.listdir("C:\\Users\\simon\\Documents\\01School\\02SOC\\SOC")
    #list_dir = os.listdir("/home/shanak/Documents/[01] Studium/SOČ/")
    for i, dir in enumerate(list_dir):
        print(f"{i}: {dir}")

    folder_number = int(input("Which folder/file you want to use? "))

    if list_dir[folder_number].endswith(".csv"):
        file_path = f"C:\\Users\\simon\\Documents\\01School\\02SOC\\SOC\\{list_dir[folder_number]}"
        #file_path = f"/home/shanak/Documents/[01] Studium/SOČ/{list_dir[folder_number]}"
        df = pd.read_csv(file_path, comment="#")
        return df
    else:
        file_path = f"C:\\Users\\simon\\Documents\\01School\\02SOC\\SOC\\{list_dir[folder_number]}"
        #file_path = f"/home/shanak/Documents/[01] Studium/SOČ/{list_dir[folder_number]}"

        file_names = os.listdir(file_path)
        file_names = [name for name in file_names if name.endswith(".csv")]
        df = pd.read_csv(f"{file_path}\\{file_names[1]}", comment="#", sep=r'\s*,\s*', engine='python')
        #df = pd.read_csv(f"{file_path}/{file_names[1]}", comment="#", sep=r'\s*,\s*', engine='python')
        for name in file_names[2:]:
            df = pd.concat([df, pd.read_csv(f"{file_path}\\{name}", comment="#", sep=r'\s*,\s*', engine='python').iloc[1:]], ignore_index=True)
            #df = pd.concat([df, pd.read_csv(f"{file_path}/{name}", comment="#", sep=r'\s*,\s*', engine='python').iloc[1:]], ignore_index=True)
        return df

df = load_data()
print(df.head())

def compute_theta(r0, ur0, r, ur, rc):
    try:
        angle = np.arctan2(ur, r - rc) - np.arctan2(ur0, r0 - rc)
    except:
        angle = 0
    return angle if angle >= 0 else angle + 2 * np.pi

def compute_center(r_values):
    r_min = np.min(r_values)
    r_max = np.max(r_values)
    return (r_min + r_max) / 2

center = compute_center(df.loc[df["init_r"] == df["init_r"].max(), "r"])
rotation_numbers = []
print("Center r =", center)
for init_r in df["init_r"].unique():
    tmp_values = df.loc[df["init_r"] == init_r, ["r", "ur"]].to_numpy()
    angle = 0
    for i in range(1, len(tmp_values)):
        angle += compute_theta(tmp_values[i-1,0], tmp_values[i-1,1], tmp_values[i,0], tmp_values[i,1], center)
    if len(tmp_values) > 1:
        rotation_numbers.append(angle / (2 * np.pi * (len(tmp_values)-1)))
    else:
        rotation_numbers.append(0)

print("Number of rotation numbers computed:", len(rotation_numbers))

print("Plotting", len(df), "points")

# fig = plt.figure(figsize=(10, 6))
# ax1 = fig.add_subplot(121)
# ax1.set_xlabel('r')
# ax1.set_ylabel('Rotation number')
# ax1.tick_params(axis='x', labelrotation=45)

# ax2 = fig.add_subplot(122)
# ax2.set_xlabel('r')
# ax2.set_ylabel('ur')
# ax2.set_title('Poincare map') 

# ax1.scatter(df['init_r'].unique(), rotation_numbers, s=4, marker='o')
# ax2.scatter(df.loc[:, "r"], df.loc[:, "ur"], s=0.1, cmap='viridis', marker='o') # c=poincare_map[:, 0]
# plt.show()

source1 = ColumnDataSource(data=dict(
    init_r=df['init_r'].unique(),
    rotation_numbers=rotation_numbers
))

source2 = ColumnDataSource(data=dict(
    r=df["r"],
    ur=df["ur"]
))

# First plot: Rotation number vs r
p1 = figure(title="Rotation numbers", x_axis_label='r', y_axis_label='Rotation number', sizing_mode="stretch_both")
p1.xaxis.major_label_orientation = 45  # rotate x-axis labels
p1.scatter(x='init_r', y='rotation_numbers',
           source=source1, size=4, marker='circle', color='navy', alpha=0.6)

# Second plot: ur vs r (Poincaré map)
p2 = figure(title="Poincaré map", x_axis_label='r', y_axis_label='ur', sizing_mode="stretch_both")
p2.scatter(x='r', y='ur', source=source2, size=3, color='green', alpha=0.7)

# Arrange the two plots side by side
layout = row(p1, p2, sizing_mode="stretch_both")
p1.add_tools(ZoomInTool(), ZoomOutTool())
p2.add_tools(ZoomInTool(), ZoomOutTool())

# Display
show(layout)

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
