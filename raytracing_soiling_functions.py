# Function Script for Raytracing-w-Soiling

import numpy as np
import ast

# From the relative positions of the receiver to the panel, the angle known as theta aim can be calculated. 
# See Aiming Strategy for LFRs for more information.
def calculate_theta_aim(Xaim, Zaim, X0, Z0):
    V = np.array([Xaim-X0, 0, Zaim - Z0])
    if V[0] == 0: Xaim += 0.00001
    theta_aim = ((Xaim-X0)/abs(Xaim-X0)) * np.arccos((Zaim - Z0) / (np.linalg.norm(V)))
    return theta_aim

def calculate_tilt(theta_transversal, theta_aim):
    return (theta_transversal + theta_aim)/2

# Returns the normal vector of a panel with a certain tilt
# This is calculated from the assumption that it is a rotation about the y axis, and so it is 
# the matrix product of the normal vector [0,0,1] and the rotation matrix Ry.
def calculate_panel_normal(tilt_rad): 
    return np.array([np.sin(tilt_rad), 0, np.cos(tilt_rad)])

# This finds the normal vector between two points
def find_normal(input_point, output_point):
    
    normal = np.array([output_point[0] - input_point[0],
                       output_point[1] - input_point[1],
                       output_point[2] - input_point[2]])
    return normal

# Distributes the parameters from the imported parameter dataframe.
def import_simulation_parameters(dataframe):
    num = dataframe.to_numpy()
    lat = float(num[0][1])
    lon = float(num[1][1])
    timezoneOffset = float(num[2][1])
    receiver_height = float(num[3][1])
    receiver_length = float(num[4][1])
    receiver_diameter = float(num[5][1])
    panel_length = float(num[6][1])
    panel_width = float(num[7][1])
    panel_height = float(num[8][1])
    panel_spacing = float(num[9][1])
    number_of_modules = int(num[10][1])
    panels_per_module = int(num[11][1])
    slope_error = float(num[12][1])
    specularity_error = float(num[13][1])

    return lat, lon, timezoneOffset, receiver_height, receiver_length, receiver_diameter, panel_length, panel_width, panel_height, panel_spacing, number_of_modules, panels_per_module, slope_error, specularity_error

