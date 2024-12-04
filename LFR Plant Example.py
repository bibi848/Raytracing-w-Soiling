"""
This script attempts to visually convey a LFR plant at a given date and time.
It models 11 LFR panels, lined up along the X (West to East) Axis, with a variable sun position, which is then plotted in an interactive 3D panel. 
XYZ: X, East. Y, North. Z, Up.
Further information for this script can be found on Aiming Strategy for Linear Fresnel Reflectors Document, and the SolTrace API documentation.
"""

# FUNCTIONS
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))

import numpy as np
import pysolar.solar as solar
from pysoltrace import PySolTrace, Point
import datetime as dt
from datetime import datetime

from raytracing_soiling_functions import calculate_theta_aim
from raytracing_soiling_functions import calculate_tilt
from raytracing_soiling_functions import calculate_panel_normal
from raytracing_soiling_functions import find_normal

from optical_geometrical_setup import op_fictitious_surface
from optical_geometrical_setup import op_cover_surface
from optical_geometrical_setup import op_heliostat_surface
from optical_geometrical_setup import op_receiver_surface
from optical_geometrical_setup import op_secondaryReflector_surface
from optical_geometrical_setup import trapezoidal_secondary_reflector

# Date and Location
# Location is Woomera and the date is set as 01/04/2023 at 11:00. This can be changed to any date.
lat, lon = -31.2,136.816667
timezoneOffset = dt.timedelta(hours = 9.5)
date = datetime(2023, 4, 1, hour=11, minute=0, second=0, tzinfo=dt.timezone(timezoneOffset))

# Plant layout
receiver_height = 4.5    # [m]
receiver_length = 12     # [m]
receiver_diameter = 0.15 # [m]
receiver_position = [0, 0, receiver_height]

panel_length = 12        # [m]
panel_width = 0.5        # [m]
panel_height = 0.5       # [m]
panel_spacing = 0.2      # [m]
# Panel Positions describes the x coordinate for each of the mirrors, ranging from -3.5 [m] to 3.75 [m],
# with a panel spacing as chosen above. 
panel_positions = np.arange(-3.5, 3.75, panel_width + panel_spacing) 
slope_error = 0.1        # [mrad]
specularity_error = 0.1  # [mrad]

stg0_height = receiver_height + 1
stg0_length = panel_length
stg0_width = panel_width * len(panel_positions) + panel_spacing * (len(panel_positions) - 1)

# Create API class instance
PT = PySolTrace()

# Finding the position of the sun
# All angle conventions are explained in Aiming Strategy for LFRs document. 
elevation_deg = solar.get_altitude(lat,lon,date) 
azimuth_deg = solar.get_azimuth(lat,lon,date) 
zenith_deg = 90 - elevation_deg  
elevation_rad, azimuth_rad, zenith_rad = (np.deg2rad(x) for x in [elevation_deg, azimuth_deg, zenith_deg])

# Describing the sun's position in terms of the azimuth and zenith. The full breakdown for this result is shown in
# Aiming Strategy for LFRs document. 
sun_position = np.array([np.sin(azimuth_rad)*np.sin(zenith_rad), np.cos(azimuth_rad)*np.sin(zenith_rad), np.cos(zenith_rad)])

# The XYZ position of the sun is then inputted into the SolTrace API. 
sun = PT.add_sun()
sun.position.x = 1000*sun_position[0]
sun.position.y = 1000*sun_position[1]
sun.position.z = 1000*sun_position[2]

# The sun position vector is normalised, to then find theta transversal and theta longitudinal. 
sn = sun_position[0:3]/np.linalg.norm(sun_position[0:3])
theta_T = np.arctan(sn[0]/sn[2])
theta_L = np.arctan(sn[1]/sn[2])

"""
Order of Stages:
0, Fictitious Surface: This is so that the sun's rays first pass through a surface before interacting with anything else, ensures the shading effect
                       of the receiver and secondary reflector are taken into account.
1, Cover: This stage contains the 'casing' for the receiver and secondary reflector. It ensures that its shading effect is taken into account.
2, Heliostats: After being shaded by the Cover stage, the light will hit the heliostats. These are angled to optimise the reflectance of the sun into the receiver.
3, Receiver & Secondary Reflector: Any rays that miss the receiver will strike the secondary reflector to hopefully increase the incident radiation on the receiver.
"""

# Stage 0, Fictitious Surface
stg0 = PT.add_stage()
stg0.is_multihit = True
stg0.is_virtual = False
stg0.name = 'Stage 0: Fictitious Surface'
stg0.position = Point(0,0,0)

optics_fictitious = op_fictitious_surface(PT, slope_error, specularity_error)
el0 = stg0.add_element()
el0.position = Point(*(7*sun_position))
aim = sun_position + 1000*find_normal(sun_position, [0,0,0])
el0.aim = Point(*aim)
el0.surface_flat()
el0.aperture_rectangle(14,14)
el0.optic = optics_fictitious
el0.interaction = 1

# Stage 1, Cover
stg1 = PT.add_stage()
stg1.is_multihit = True
stg1.is_tracethrough = True
stg1.name = 'Stage 1: Cover'
stg1.position = Point(0,0,0)
optics_cover = op_cover_surface(PT, slope_error, specularity_error)
el1 = trapezoidal_secondary_reflector(stg1, optics_cover, receiver_height, receiver_length)

# Setting up the Heliostats
stg2 = PT.add_stage()
stg2.is_multihit = True
stg2.name = 'Stage 2: Heliostats'
stg2.position = Point(0,0,0)
for p in range(len(panel_positions)):
    optics_heliostat_p = op_heliostat_surface(PT, slope_error, specularity_error, 1, p)

    heliostat_position = [panel_positions[p], 0, panel_height]
    # From the relative positions of the panel to the receiver, the panel's tilt and normal are calculated.
    theta_aim = calculate_theta_aim(Xaim=receiver_position[0], Zaim=receiver_position[2], X0=heliostat_position[0], Z0=heliostat_position[2])
    tilt = calculate_tilt(theta_T, theta_aim)
    panel_normal = calculate_panel_normal(tilt)

    el2 = stg2.add_element()
    el2.position = Point(*heliostat_position)
    aim = heliostat_position + 1000*panel_normal
    el2.aim = Point(*aim)
    el2.optic = optics_heliostat_p
    el2.surface_flat()
    el2.aperture_rectangle(panel_width, panel_length)

# Receiver & Secondary Reflector
stg3 = PT.add_stage()
stg3.is_multihit = True
stg3.name = 'Stage 3: Receiver & Secondary Reflector'
stg3.position = Point(0,0,0)
optics_receiver = op_receiver_surface(PT)

el3 = stg3.add_element()
el3.position = Point(*receiver_position)
el3.aim = Point(0,0,0)
el3.optic = optics_receiver
el3.surface_cylindrical(receiver_diameter/2)
el3.aperture_singleax_curve(0, 0, receiver_length) # (inner coordinate of revolved section, outer coordinate of revolved section, 
                                                #  length of revolved section along axis of revolution)
optics_secondary = op_secondaryReflector_surface(PT, slope_error, specularity_error)
el3 = trapezoidal_secondary_reflector(stg3, optics_secondary, receiver_height, receiver_length)


# Simulation Parameters
PT.num_ray_hits = 1e5
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True
PT.is_surface_errors = True

if __name__ == "__main__":

    PT.run(-1, True, 8)
    print("Num rays traced: {:d}".format(PT.raydata.index.size))
    PT.plot_trace()