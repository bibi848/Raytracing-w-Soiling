# This script will graph two sets of data together, taking in the file path and header names as inputs.

import pandas as pd
import matplotlib.pyplot as plt
import os

# Replace the folder_path and file_name with those that you need.
folder_path = 'C:/Users/oscar/Documents/VRES/Raytracing-w-Soiling/'
file_name = 'COMMERCIAL_lookup_table-SOILED-14days.csv'
header_x = 'rho_avg_14'
header_y = 'dni'

file_path = os.path.join(folder_path, file_name)

def plot_from_csv(file_path, header_x, header_y, marker = 'o', linestyle = '-', color = 'b'):
    data = pd.read_csv(file_path)
    x_data = data[header_x]
    y_data = data[header_y]

    plt.figure(figsize=(8, 6))
    plt.plot(x_data, y_data, marker, linestyle, color)
    plt.title(f"{header_y} vs {header_x}")
    plt.xlabel(header_x)
    plt.ylabel(header_y)
    plt.grid(True)
    plt.tight_layout()
    plt.show()