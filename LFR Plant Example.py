# This script attempts to visually convey a LFR plant at a given date and time.
# It models 11 LFR panels, lined up along the X (West to East) Axis, with a variable sun position, which is then 
# plotted in an interactive 3D panel. 
# XYZ: X, East. Y, North. Z, Up.
# Further information for this script can be found on Aiming Strategy for Linear Fresnel Reflectors Document, and the 
# SolTrace API documentation.

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

# Setting up the Heliostats
st1 = PT.add_stage()
st1.position = Point(0,0,0)
for p in range(len(panel_positions)):
    optics_heliostat_p = PT.add_optic("Heliostat_{p}")
    optics_heliostat_p.front.reflectivity = 1.0

    heliostat_position = [panel_positions[p], 0, panel_height]
    # From the relative positions of the panel to the receiver, the panel's tilt and normal are calculated.
    theta_aim = calculate_theta_aim(Xaim=receiver_position[0], Zaim=receiver_position[2], X0=heliostat_position[0], Z0=heliostat_position[2])
    tilt = calculate_tilt(theta_T, theta_aim)
    panel_normal = calculate_panel_normal(tilt)

    el = st1.add_element()
    el.position = Point(*heliostat_position)
    aim = heliostat_position + 1000*panel_normal
    el.aim = Point(*aim)
    el.optic = optics_heliostat_p
    el.surface_flat()
    el.aperture_rectangle(panel_width, panel_length)

# Receiver
optics_receiver = PT.add_optic("Receiver")
optics_receiver.front.reflectivity = 0.0

sta = PT.add_stage()
ela = sta.add_element()
ela.position = Point(*receiver_position)
ela.aim = Point(0,0,0)
ela.optic = optics_receiver
ela.surface_cylindrical(receiver_diameter/2)
ela.aperture_singleax_curve(0, 0, receiver_length) # (inner coordinate of revolved section, outer coordinate of revolved section, 
                                                   #  length of revolved section along axis of revolution)


# Simulation Parameters
PT.num_ray_hits = 1e6
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True
PT.is_surface_errors = True

if __name__ == "__main__":

    PT.run(-1, True, 8)
    print("Num rays traced: {:d}".format(PT.raydata.index.size))
    PT.plot_trace()