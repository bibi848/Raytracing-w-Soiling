
# Function Script for Raytracing-w-Soiling

import numpy as np
from scipy.spatial import ConvexHull

# From the relative positions of the receiver to the panel, the angle known as theta aim can be calculated. 
# See Aiming Strategy for LFRs for more information.
def calculate_theta_aim(Xaim, Zaim, X0, Z0):
    V = np.array([Xaim-X0, 0, Zaim - Z0])
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

