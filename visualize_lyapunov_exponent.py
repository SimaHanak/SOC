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
        #df = pd.read_csv(f"{file_path}/{file_names[0]}", comment="#", sep=r'\s*,\s*', engine='python')
        return df[1:]

df = load_data()
print(df.head())

source1 = ColumnDataSource(data=dict(
    logt=np.log10(df["t"]),
    lyapunov=np.log10(np.abs(df["sum_log_stretch"] / df["t"]))
))

p1 = figure(title="Lyapunov exponent", x_axis_label='log(t)', y_axis_label='log(Lyapunov exponent)', sizing_mode="stretch_both")
p1.axis.axis_label_text_font_size = "14pt"
p1.axis.major_label_text_font_size = "12pt"
p1.title.text_font_size = "16pt"
p1.scatter(x='logt', y='lyapunov',
           source=source1, size=1, marker='circle', alpha=0.3)

p1.add_tools(ZoomInTool(), ZoomOutTool())

# Display
show(p1)