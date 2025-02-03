# %%

import matplotlib.pyplot as plt
from numpy import pi,deg2rad,abs,cos,sin,sqrt,linspace,zeros,arccos,r_,argmin
from scipy.optimize import bisect

""" Used equations from [1] to obtain the ideal profile for the CPC.
[1] A. Rabl, “Solar concentrators with maximal concentration for cylindrical absorbers,” Appl. Opt., AO, vol. 15, no. 7, pp. 1871–1873, Jul. 1976, doi: 10.1364/AO.15.001871.
"""

a = 0.07/2.0 # absorber radius [m]
a_design = a
cpc_depth = 0.11
theta_a = deg2rad(56)

def aperature_angle(aw, d):

    if d < a:
        raise ValueError("Depth of collector must be at least on absorber radius.")

    xf,yf = aw,d
    fun = lambda t: xf*cos(t)+yf*sin(t)-a
    t = bisect(fun,0,pi)
    xc,yc = a*cos(t),a*sin(t)
    dp = ( (xf-xc)*0 + (yf-yc)*d )/sqrt(xf**2+yf**2)/sqrt(d**2)
    theta_a = arccos(dp)

    return theta_a,t

def distance_from_tangent(θ,theta_a,ad):
    if abs(θ) <= theta_a + pi/2:
        ρ = ad*θ

    elif (abs(θ) > theta_a + pi/2) and (abs(θ) <= 3.0*pi/2-theta_a):
        ρ = ad* (θ+theta_a+pi/2 - cos(θ-theta_a))/(1+sin(θ-theta_a))

    else:
        raise ValueError("\theta too large.")

    return ρ

def reflector_point(θ,ρ,ad):

    # find point on absorber surface & compute tangent unit vector
    xa,ya = ad*cos(θ-pi/2),ad*sin(θ-pi/2)
    L = sqrt(xa**2+ya**2)
    xt,yt = ya/L,-xa/L

    return xa+ρ*xt,ya+ρ*yt

def design_secondary_concentrator(theta_a,theta_max,ad,n_thetas=100,clip_points=0):
    theta = linspace(0,theta_max,n_thetas)
    x,y = zeros(theta.shape),zeros(theta.shape)

    for ii,t in enumerate(theta):
        r = distance_from_tangent(t,theta_a,ad)
        x[ii],y[ii] = reflector_point(t,r,ad)

    x = x[clip_points::]
    y = y[clip_points::]
    x = r_[-x[::-1],x]
    y = r_[y[::-1],y]

    return x,y

# Sampling the CPC curve
# in: x and y coordinates of curve and number of divisions, out: midpoints of panels, start points, end points

def find_midpoint(x1, x2, y1, y2):
    x_m = (x1+x2)/2
    y_m = (y1+y2)/2
    return x_m, y_m


def sample_CPC(x, y, number_divisions):

    x_shift = int(len(x)/number_divisions)
    i = 0
    sample_points_x = []
    sample_points_y = []
    while i < len(x)/2:
        sample_points_x.append(x[i])
        sample_points_y.append(y[i])
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

a = 0.15/2.0 # absorber radius [m]
a_design = a
cpc_depth = 0.11
theta_a = deg2rad(70)

x, y = design_secondary_concentrator(theta_a, 3*pi/2-theta_a, a_design, n_thetas=100, clip_points=20)

sample_points_x, sample_points_y, midpoints_x, midpoints_y = sample_CPC(x, -1*y, 12)

fig,ax = plt.subplots()

ax.plot(a*cos(linspace(0,2*pi,1000)),a*sin(linspace(0,2*pi,1000)),'b-',label='receiver tube')
ax.plot(x,-y,'k-',label='secondary collector')
ax.scatter(sample_points_x, sample_points_y, color = 'red', label = 'samples')
ax.plot(sample_points_x, sample_points_y, color = 'green')
ax.scatter(midpoints_x, midpoints_y, color = 'green', label = 'midpoints')

ax.axis('equal')
ax.legend()
# %%
