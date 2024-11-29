
# This script attempts to visually convey a LFR plant with a trapezoidal secondary reflector at a given date and time.

# System Setup
# Coordinate System, XYZ: X = West, Y = Perpendicular to the ground, Z = North.
# The 1st mirror is the East-most mirror, and the 11th mirror is the West-most mirror.

# Load the pysoltrace api from the parent directory ---
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))
# ----------

from pysoltrace import PySolTrace, Point
from math import sin,cos, pi
import random



# Create API class instance
PT = PySolTrace()

# Plant layout
receiver_height = 4.5    # [m]
receiver_length = 12.18  # [m]
receiver_diameter = 0.15 # [m]
receiver_position = Point(0.0, 0.0, receiver_height)

panel_length = 12.18 # [m]
panel_width = 0.5    # [m]

# Optics
opt_ref = PT.add_optic("Reflector") 
opt_ref.front.reflectivity = 1.0

opt_abs = PT.add_optic("Absorber")
opt_abs.front.reflectivity = 0.0

# Sun
sun = PT.add_sun()
sun.position.x = 0.0
sun.position.y = 0.0
sun.position.z = 100.0

receiver_pos = Point(0.0, 0.0, 10)

# Reflector
st = PT.add_stage()
el = st.add_element()
el.position.x = 2.0
el.position.y = 1.5

r_vector = (receiver_position - el.position).unitize()
s_vector = sun.position.unitize()
avg_vector = (r_vector + s_vector)/2

el.aim = el.position + avg_vector*100
el.optic = opt_ref
el.zrot = PT.util_calc_zrot_azel(avg_vector)
el.surface_flat()
el.aperture_rectangle(panel_length, panel_width)

# Absorber
sta = PT.add_stage()
ela = sta.add_element()
ela.position = receiver_position
ela.aim = Point(0,100,0)
ela.optic = opt_abs
ela.surface_cylindrical(receiver_diameter/2)
ela.aperture_singleax_curve(0, 0, receiver_length)

# Simulation Parameters
PT.num_ray_hits = 1e6  
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True
PT.is_surface_errors = True

if __name__ == "__main__":

    PT.run(-1, True, 8)
    print("Num rays traced: {:d}".format(PT.raydata.index.size))
    PT.plot_trace()
