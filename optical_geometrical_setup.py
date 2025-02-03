# Optical and Geometrical Setup for PySolTrace

from pysoltrace import Point
import numpy as np


# Optical
def stg_surface(PT, name, f_reflectivity, b_reflectivity, f_transmissivity, b_transmissivity, slope_error = 0, specularity_error = 0):
    optics = PT.add_optic(name)
    optics.front.reflectivity = f_reflectivity
    optics.back.reflectivity = b_reflectivity
    optics.front.transmissivity = f_transmissivity
    optics.back.transmissivity = b_transmissivity
    optics.front.slope_error = slope_error
    optics.front.specularity_error = specularity_error
    optics.back.slope_error = slope_error
    optics.back.specularity_error = specularity_error
    return optics

def op_fictitious_surface(PT, slope_error, specularity_error):
    optics_fictitious = PT.add_optic("Fictitious Surface")
    optics_fictitious.front.reflectivity = 0.0 
    optics_fictitious.back.reflectivity = 0.0
    optics_fictitious.front.transmissivity = 1.0
    optics_fictitious.back.transmissivity = 1.0
    optics_fictitious.front.refraction_real = 1.0
    optics_fictitious.back.refraction_real = 1.0 
    optics_fictitious.front.slope_error = slope_error
    optics_fictitious.front.specularity_error = specularity_error
    optics_fictitious.back.slope_error = slope_error
    optics_fictitious.back.specularity_error = specularity_error
    return optics_fictitious

def op_cover_surface(PT, slope_error, specularity_error):
    optics_secondary = PT.add_optic("Cover")
    optics_secondary.front.reflectivity = 1.0 
    optics_secondary.front.transmissivity = 0.0
    optics_secondary.back.reflectivity = 1.0 
    optics_secondary.back.transmissivity = 0.0
    optics_secondary.front.slope_error = slope_error
    optics_secondary.front.specularity_error = specularity_error
    optics_secondary.back.slope_error = slope_error
    optics_secondary.back.specularity_error = specularity_error
    return optics_secondary

def op_heliostat_surface(PT, slope_error, specularity_error, reflectivity, p):
    optics_p = PT.add_optic(f"Heliostat_{p}") 
    optics_p.front.reflectivity = reflectivity  
    optics_p.front.transmissivity = 0.0
    optics_p.front.slope_error = slope_error
    optics_p.front.specularity_error = specularity_error
    optics_p.back.reflectivity = 0.0
    optics_p.back.transmissivity = 0.0
    return optics_p

def op_receiver_surface(PT):
    optics_receiver = PT.add_optic("Receiver")
    optics_receiver.front.reflectivity = 0.06
    optics_receiver.back.reflectivity = 0.06
    optics_receiver.front.transmissivity = 0.0
    optics_receiver.back.transmissivity = 0.0
    return optics_receiver

def op_secondaryReflector_surface(PT, slope_error, specularity_error):
    optics_secondary = PT.add_optic("Secondary Reflector")
    optics_secondary.front.reflectivity = 0.9
    optics_secondary.front.transmissivity = 0.0
    optics_secondary.back.reflectivity = 0.9
    optics_secondary.back.transmissivity = 0.0
    optics_secondary.front.slope_error = slope_error
    optics_secondary.front.specularity_error = specularity_error
    optics_secondary.back.slope_error = slope_error
    optics_secondary.back.specularity_error = specularity_error
    return optics_secondary

def trapezoidal_secondary_reflector(stg, optics, receiver_height, receiver_length):

    spos1 = [-0.15, 0, receiver_height-0.06]
    spos2 = [0.15, 0, receiver_height-0.06]
    spos3 = [0, 0, receiver_height+0.04]

    el = stg.add_element() 
    el.optic = optics
    el.position = Point(*spos1) 
    el.aim = Point(-300, 0, 275) 
    el.aperture_rectangle(0.27, receiver_length); 
    el.surface_flat(); 
    el.interaction = 2 

    el = stg.add_element() 
    el.optic = optics
    el.position = Point(*spos2) 
    el.aim = Point(300, 0, 275) 
    el.aperture_rectangle(0.27, receiver_length)
    el.surface_flat(); 
    el.interaction = 2 

    el = stg.add_element() 
    el.optic = optics
    el.position = Point(*spos3)  
    el.aim = Point(0.000001,0,0)
    el.aperture_rectangle(0.12, receiver_length)
    el.surface_flat()
    el.interaction = 2

    return el


# Design of a CPC Reflector

""" Used equations from [1] to obtain the ideal profile for the CPC.
[1] A. Rabl, “Solar concentrators with maximal concentration for cylindrical absorbers,” Appl. Opt., AO, vol. 15, no. 7, pp. 1871–1873, Jul. 1976, doi: 10.1364/AO.15.001871.

Functions: aperature_angle(...), distance_from_tangent(...), reflector_point(...), design_secondary_reflector(...) were provided by Prof. Michael Cholette.
"""

from scipy.optimize import bisect

def aperature_angle(aw, d, a):

    if d < a:
        raise ValueError("Depth of collector must be at least on absorber radius.")

    xf, yf = aw, d
    fun = lambda t: xf*np.cos(t) + yf*np.sin(t) - a
    t = bisect(fun, 0, np.pi)
    xc, yc = a*np.cos(t), a*np.sin(t)
    dp = ((xf-xc) * 0 + (yf - yc) *d) / np.sqrt(xf**2 + yf**2) / np.sqrt(d**2)
    theta_a = np.arccos(dp)

    return theta_a,t

def distance_from_tangent(θ, theta_a, ad):
    if abs(θ) <= theta_a + np.pi/2:
        rho = ad*θ

    elif (abs(θ) > theta_a + np.pi/2) and (abs(θ) <= 3.0*np.pi/2-theta_a):
        rho = ad* (θ+theta_a + np.pi/2 - np.cos(θ-theta_a)) / (1+np.sin(θ-theta_a))
    else:
        raise ValueError("\theta too large.")

    return rho

def reflector_point(θ, rho, ad):

    xa,ya = ad*np.cos(θ - np.pi/2), ad*np.sin(θ - np.pi/2)
    L = np.sqrt(xa**2+ya**2)
    xt, yt = ya/L, -xa/L

    return xa+rho*xt, ya+rho*yt

def design_secondary_concentrator(theta_a, theta_max,ad, n_thetas=100, clip_points=0):
    theta = np.linspace(0, theta_max, n_thetas)
    x,y = np.zeros(theta.shape), np.zeros(theta.shape)

    for ii,t in enumerate(theta):
        r = distance_from_tangent(t,theta_a,ad)
        x[ii],y[ii] = reflector_point(t,r,ad)

    x = x[clip_points::]
    y = y[clip_points::]
    x = np.r_[-x[::-1],x]
    y = np.r_[y[::-1],y]

    return x,y

# Finds the midpoint of two points on the same plane
def find_midpoint(x1, x2, y1, y2):
    x_m = (x1+x2)/2
    y_m = (y1+y2)/2
    return x_m, y_m

# From the array of coordinates making up the CPC curve, the curve is discretized to be used in 
# the raytracing.
def sample_CPC(x, y, number_divisions, height_offset):

    x_shift = int(len(x)/number_divisions)
    i = 0
    sample_points_x = []
    sample_points_y = []
    while i < len(x)/2:
        sample_points_x.append(x[i])
        sample_points_y.append(y[i] + height_offset)
        i += x_shift

    reverse_y = sample_points_y[::-1]
    reverse_x = sample_points_x[::-1]

    for i in range(len(reverse_x)):

        sample_points_x.append(-1*reverse_x[i])
        sample_points_y.append(reverse_y[i])
    
    midpoints_x = []
    midpoints_y = []
    for i in range(len(sample_points_x) - 1):
        x_m, y_m = find_midpoint(sample_points_x[i], sample_points_x[i+1], sample_points_y[i], sample_points_y[i+1])
        midpoints_x.append(x_m)
        midpoints_y.append(y_m)

    return sample_points_x, sample_points_y, midpoints_x, midpoints_y

def CPC_positioning(stg, optics, panel_length, sample_points_x, sample_points_y, midpoints_x, midpoints_y):

    for p in range(len(midpoints_x)):

        reflector_position = np.array([midpoints_x[p], 0, midpoints_y[p]])

        point1 = reflector_position
        point2 = np.array([sample_points_x[p], 6, sample_points_y[p]])
        point3 = np.array([sample_points_x[p+1], 6, sample_points_y[p+1]])

        v1 = point2 - point1
        v2 = point3 - point1

        n = np.cross(v1, v2)
        width = np.linalg.norm(point3-point2)

        el = stg.add_element()
        el.position = Point(*reflector_position)
        aim = reflector_position + 1000*n
        el.aim = Point(*aim)
        el.optic = optics
        el.surface_flat()
        el.aperture_rectangle(width, panel_length)