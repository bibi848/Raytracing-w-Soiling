
# This script attempts to visually convey a LFR plant with a trapezoidal secondary reflector at a given date and time.

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
# Location is Woomera and the date is set as 01/04/2023 at 09:00
lat, lon = -31.2,136.816667
timezoneOffset = dt.timedelta(hours = 9.5)
date = datetime(2023, 4, 1, hour=9, minute=0, second=0, tzinfo=dt.timezone(timezoneOffset))

# Create API class instance
PT = PySolTrace()

# Find the position of the sun
elevation_deg = solar.get_altitude(lat,lon,date) 
azimuth_deg = solar.get_azimuth(lat,lon,date) 
zenith_deg = 90 - elevation_deg  
elevation_rad, azimuth_rad, zenith_rad = (np.deg2rad(x) for x in [elevation_deg, azimuth_deg, zenith_deg])

sun_position = np.array([-np.sin(azimuth_rad)*np.sin(zenith_rad), np.cos(zenith_rad), np.cos(azimuth_rad)*np.sin(zenith_rad)])

sun = PT.add_sun()
sun.position.x = 1000*sun_position[0] 
sun.position.y = 1000*sun_position[1]
sun.position.z = 1000*sun_position[2]


# Plant layout
receiver_height = 4.5    # [m]
receiver_length = 12.18  # [m]
receiver_diameter = 0.15 # [m]
receiver_position = [0, 0, receiver_height]

panel_length = 12.18 # [m]
panel_width = 0.5    # [m]
panel_height = 0.5   # [m]
heliostat_position = [2, 0, panel_height]

# Tilt of panel
sn = sun_position[0:2]/np.linalg.norm(sun_position[0:2])
theta_T = np.arctan(sn[0]/sn[1])
theta_L = np.arctan(np.cos(azimuth_rad) * np.tan(zenith_rad))
theta_aim = calculate_theta_aim(0, receiver_height, 2, panel_height) # Change to receiver position [0]
tilt = calculate_tilt(theta_T, theta_aim)
panel_normal = calculate_panel_normal(tilt)

# Optical & Geometrical Setup
# Heliostat
optics_heliostat = PT.add_optic("Heliostat")
optics_heliostat.front.reflectivity = 1.0

st1 = PT.add_stage()
st1.position = Point(0,0,0)
el = st1.add_element()
el.position = Point(*heliostat_position)

aim = heliostat_position + 1000*panel_normal
el.aim = Point(*aim)

el.optic = optics_heliostat
el.surface_flat()
el.aperture_rectangle(panel_length, panel_width)

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

#print(PT.util_calc_zrot_azel(avg_vector))

if __name__ == "__main__":

    PT.run(-1, True, 8)
    print("Num rays traced: {:d}".format(PT.raydata.index.size))
    PT.plot_trace()

