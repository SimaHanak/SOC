import pandas as pd
import numpy as np
import os
from bokeh.plotting import figure, show, ColumnDataSource
from bokeh.layouts import row
from bokeh.models import ZoomInTool, ZoomOutTool
from math import sin, cos

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
        df = pd.read_csv(f"{file_path}\\{file_names[0]}", comment="#", sep=r'\s*,\s*', engine='python')
        #df = pd.read_csv(f"{file_path}/{file_names[1]}", comment="#", sep=r'\s*,\s*', engine='python')
        return df[:100000]

df = load_data()
print(df.head())

df_potential = pd.read_csv("C:\\Users\\simon\\Documents\\01School\\02SOC\\SOC\\V_eff_0.csv", comment="#")

source1 = ColumnDataSource(data=dict(
    r=pd.concat([df['r'], df_potential['r']]),
    z=pd.concat([df['z'], df_potential['z']]),
    colors=['navy'] * len(df) + ['red'] * len(df_potential),
))
source2 = ColumnDataSource(data=dict(
    rsin=[r*sin(phi) for r, phi in zip(df['r'], df['phi'])],
    rcos=[r*cos(phi) for r, phi in zip(df['r'], df['phi'])],
))

# First plot: Rotation number vs r
p1 = figure(title="Chaotická trajektorie v rovině rz", x_axis_label='r', y_axis_label='z', sizing_mode="stretch_both")
p1.axis.axis_label_text_font_size = "14pt"
p1.axis.major_label_text_font_size = "12pt"
p1.title.text_font_size = "16pt"
p1.scatter(x='r', y='z',
           source=source1, size=1, marker='circle', color='colors', alpha=0.3)

# Second plot: ur vs r (Poincaré map)
p2 = figure(title="Chaotická trajektorie v ekvatoriální rovině", x_axis_label='rsin(φ)', y_axis_label='rcos(φ)', sizing_mode="stretch_both")
p2.axis.axis_label_text_font_size = "14pt"
p2.axis.major_label_text_font_size = "12pt"
p2.title.text_font_size = "16pt"
p2.scatter(x='rsin', y='rcos', source=source2, size=1, color='green', alpha=0.3)

# Arrange the two plots side by side
layout = row(p1, p2, sizing_mode="stretch_both")
p1.add_tools(ZoomInTool(), ZoomOutTool())
p2.add_tools(ZoomInTool(), ZoomOutTool())

# Display
show(layout)