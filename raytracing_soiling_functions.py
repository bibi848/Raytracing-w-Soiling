
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

# From the plant layount, this finds the equivalent volume as a cube
def find_equivalent_cube(plength, pwidth, pspacing, pheight, no_panels):
    twidth =  pwidth * no_panels + pspacing * ((no_panels) - 1)
    theight = pwidth/2 + pheight
    cube = np.array([
        [twidth/2, -plength/2, 0],
        [twidth/2, -plength/2, theight],
        [-twidth/2, -plength/2, 0],
        [-twidth/2, -plength/2, theight],
        [-twidth/2, plength/2, 0],
        [-twidth/2, plength/2, theight],
        [twidth/2, plength/2, 0],
        [twidth/2, plength/2, theight],
    ])
    return cube

# From the cube of the solar field, its projection points are mapped to a plane about the reference point
def project_to_plane(points, sn, ref_point):
    projection_matrix = np.eye(3) - np.outer(sn, sn)  # Projection matrix
    projected_points = [projection_matrix @ (p - ref_point) + ref_point for p in points]
    return np.array(projected_points)

def points_to_2_quadrilaterals(projected_points):

    # 0, 1, 3, 5, 6, 7
    # first quad: 0,1,3,7
    # 2nd quad: 3, 5, 6, 7
    quad1 = np.array([projected_points[0], projected_points[1], projected_points[3], projected_points[6]])
    quad2 = np.array([projected_points[3], projected_points[4], projected_points[5], projected_points[6]])

    return quad1, quad2

def fic_surf_pos(sun_position, field_coords, distance_multiplier):
    surface_coords = []
    for coord in field_coords:
        res = coord + sun_position*distance_multiplier
        surface_coords.append(res)
    return surface_coords