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
import matplotlib.pyplot as plt
import pysolar.solar as solar
from pysoltrace import PySolTrace, Point
import datetime as dt
from datetime import datetime

from raytracing_soiling_functions import calculate_theta_aim
from raytracing_soiling_functions import calculate_tilt
from raytracing_soiling_functions import calculate_panel_normal
from raytracing_soiling_functions import find_normal
from raytracing_soiling_functions import find_equivalent_cube
from raytracing_soiling_functions import project_to_plane
from raytracing_soiling_functions import points_to_2_quadrilaterals
from raytracing_soiling_functions import fic_surf_pos

from optical_geometrical_setup import op_fictitious_surface
from optical_geometrical_setup import op_cover_surface
from optical_geometrical_setup import op_heliostat_surface
from optical_geometrical_setup import op_receiver_surface
from optical_geometrical_setup import op_secondaryReflector_surface
from optical_geometrical_setup import trapezoidal_secondary_reflector

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
equivalent_cube = find_equivalent_cube(panel_length, panel_width, panel_spacing, panel_height, len(panel_positions))
slope_error = 0.1        # [mrad]
specularity_error = 0.1  # [mrad]

stg0_height = receiver_height
stg0_length = panel_length 
stg0_width = panel_width * len(panel_positions) + panel_spacing * (len(panel_positions) - 1) 
distance_multiplier = 1 # Scaling factor which pushes the fictitious surface away from the solar field
field_coords = [[-3.5, 6], [-3.5, -6], [3.75, 6], [3.75, -6]]

# Create API class instance
PT = PySolTrace()

# Describing the sun's position in terms of the azimuth and zenith. The full breakdown for this result is shown in
# Aiming Strategy for LFRs document. 
sun_position = np.array([0.01, 0.01, 10])

# The XYZ position of the sun is then inputted into the SolTrace API. 
sun = PT.add_sun()
sun.position.x = 2000*sun_position[0]
sun.position.y = 2000*sun_position[1]
sun.position.z = 2000*sun_position[2]

# The sun position vector is normalised, to then find theta transversal and theta longitudinal. 
sn = sun_position[0:3]/np.linalg.norm(sun_position[0:3])
theta_T = np.arctan(sn[0]/sn[2])
theta_L = np.arctan(sn[1]/sn[2])

# Stage 0, Fictitious Surface
# stg0 = PT.add_stage()
# stg0.is_multihit = True
# stg0.is_virtual = False
# stg0.name = 'Stage 0: Fictitious Surface'
# stg0.position = Point(0,0,0)

# optics_fictitious = op_fictitious_surface(PT, slope_error, specularity_error)
# el0 = stg0.add_element()
# el0.position = Point(distance_multiplier*sun_position[0], distance_multiplier*sun_position[1], distance_multiplier*sun_position[2])
# aim = sun_position + 1000*find_normal(sun_position, [0,0,0])
# el0.aim = Point(*aim)
# el0.surface_flat()
# el0.aperture_rectangle(stg0_width, stg0_length)
# el0.optic = optics_fictitious
# el0.interaction = 1


stg1 = PT.add_stage()
stg1.is_multihit = True
stg1.name = 'Stage 0: Fictitious Surface'
stg1.position = Point(0,0,0)

optics_heliostat_p = op_heliostat_surface(PT, slope_error, specularity_error, 1, 1)
el1 = stg1.add_element()
#el1.position = Point(0,0,0)
aim = sun_position + 1000*find_normal(sun_position, [0,0,0])
el1.aim = Point(*aim)
el1.aperture_quadrilateral(0,0,0,1,-1,0,-1,1)
# el1.aperture_rectangle(stg0_width, stg0_length)
el1.surface_flat()
el1.optic = optics_heliostat_p


# Simulation Parameters
PT.num_ray_hits = 1e4
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True
PT.is_surface_errors = True

if __name__ == "__main__":

    PT.run(-1, True, 4)
    print("Num rays traced: {:d}".format(PT.raydata.index.size))
    PT.plot_trace()