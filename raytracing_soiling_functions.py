
# Function Script for Raytracing-w-Soiling

import numpy as np

def calculate_theta_aim(Xaim, Zaim, X0, Z0):
    V = np.array([Xaim-X0, 0, Zaim - Z0])
    theta_aim = ((Xaim-X0)/abs(Xaim-X0)) * np.arccos((Zaim - Z0) / (np.linalg.norm(V)))
    return theta_aim

def calculate_tilt(theta_transversal, theta_aim):
    return (theta_transversal + theta_aim)/2

def calculate_panel_normal(tilt_rad): # Returns the normal vector of a panel with a certain tilt
    return np.array([np.sin(tilt_rad), np.cos(tilt_rad), 0])

    