"""
This script attempts to visually convey a LFR plant at a given date and time.
It models 11 LFR panels, lined up along the X (West to East) Axis, with a variable sun position, which is then plotted in an interactive 3D panel. 
XYZ: X, East. Y, North. Z, Up.
Further information for this script can be found on Aiming Strategy for Linear Fresnel Reflectors Document, and the SolTrace API documentation.
"""

# FUNCTIONS
import os
import numpy as np
import pandas as pd
import pysolar.solar as solar
from pysoltrace import PySolTrace, Point
import datetime as dt
from datetime import datetime
import pysolar.radiation as radiation
import warnings
warnings.filterwarnings("ignore", message="no explicit representation of timezones available for np.datetime64")

from raytracing_soiling_functions import calculate_theta_aim
from raytracing_soiling_functions import calculate_tilt
from raytracing_soiling_functions import calculate_panel_normal
from raytracing_soiling_functions import import_simulation_parameters

from optical_geometrical_setup import op_fictitious_surface
from optical_geometrical_setup import op_cover_surface
from optical_geometrical_setup import op_heliostat_surface
from optical_geometrical_setup import op_receiver_surface
from optical_geometrical_setup import op_secondaryReflector_surface
from optical_geometrical_setup import trapezoidal_secondary_reflector

# Importing plant parameters
current_script_path = os.path.abspath(__file__)           
current_directory = os.path.dirname(current_script_path)  
csv_path = os.path.join(current_directory, "CSV Files/Simulation Parameters.csv")
lat, lon, hour_offset, receiver_height, receiver_length, receiver_diameter, panel_length, panel_width, panel_height, panel_spacing, panels_min_max, slope_error, specularity_error = import_simulation_parameters(pd.read_csv(csv_path))

# Date and Location
# Location is Woomera and the date is set as 01/04/2018 at 11:00. This can be changed to any date.
timezoneOffset = dt.timedelta(hours = hour_offset)
date = datetime(2018, 4, 1, hour=11, minute=0, second=0, tzinfo=dt.timezone(timezoneOffset))

# Plant layout
receiver_position = [0, 0, receiver_height]

panel_positions = np.arange(panels_min_max[0], panels_min_max[1], panel_width + panel_spacing) # Describes the x-coordinate for each of the mirrors, 
                                                                                               # ranging from -3.5 [m] to 3.75 [m].
stg1_length = panel_length 
stg1_width = panel_width * len(panel_positions) + panel_spacing * (len(panel_positions) - 1) 
distance_multiplier = 10                               # Scaling factor which pushes the fictitious surface away from the solar field.
x_shift = (panel_positions[0] + panel_positions[-1])/2 # As the solar field is not exactly centered along the x-axis, there is a shift 
                                                       # required for the aiming algorithm.

# Create API class instance
PT = PySolTrace()

# Finding the position of the sun
# All angle conventions are explained in Aiming Strategy for LFRs document. 
elevation_deg = solar.get_altitude(lat,lon,date) 
azimuth_deg = solar.get_azimuth(lat,lon,date) 
zenith_deg = 90 - elevation_deg  
elevation_rad, azimuth_rad, zenith_rad = (np.deg2rad(x) for x in [elevation_deg, azimuth_deg, zenith_deg])
solar_radiation = radiation.get_radiation_direct(date,elevation_deg)

# Describing the sun's position in terms of the azimuth and zenith. The full breakdown for this result is shown in
# Aiming Strategy for LFRs document. 
sun_position = np.array([np.sin(azimuth_rad)*np.sin(zenith_rad), np.cos(azimuth_rad)*np.sin(zenith_rad), np.cos(zenith_rad)])

# The XYZ position of the sun is then inputted into the SolTrace API. 
sun = PT.add_sun()
sun.position.x = sun_position[0]
sun.position.y = sun_position[1]
sun.position.z = sun_position[2]

# The sun position vector is normalised, to then find theta transversal and theta longitudinal. 
sn = sun_position[0:3]/np.linalg.norm(sun_position[0:3])
theta_T = np.arctan(sn[0]/sn[2])
theta_L = np.arctan(sn[1]/sn[2])

"""
Order of Stages:
1, Fictitious Surface: This is so that the sun's rays first pass through a surface before interacting with anything else, ensures the shading effect
                       of the receiver and secondary reflector are taken into account.
2, Cover: This stage contains the 'casing' for the receiver and secondary reflector. It ensures that its shading effect is taken into account.
3, Heliostats: After being shaded by the Cover stage, the light will hit the heliostats. These are angled to optimise the reflectance of the sun into the receiver.
4, Receiver & Secondary Reflector: Any rays that miss the receiver will strike the secondary reflector to hopefully increase the incident radiation on the receiver.
"""

# Stage 1, Fictitious Surface
stg1 = PT.add_stage()
stg1.is_virtual = True
stg1.name = 'Stage 1: Fictitious Surface'
stg1.position = Point(0,0,0)

optics_fictitious = op_fictitious_surface(PT, slope_error, specularity_error)
el1 = stg1.add_element()
el1.position = Point(distance_multiplier*(sun_position[0]+x_shift), 
                     distance_multiplier*sun_position[1], 
                     distance_multiplier*sun_position[2])
el1.aim = Point(distance_multiplier*(sun_position[0]+x_shift), distance_multiplier*sun_position[1], 0)
el1.surface_flat()
el1.aperture_rectangle(stg1_width, stg1_length)
el1.optic = optics_fictitious

# Stage 2, Cover
stg2 = PT.add_stage()
stg2.is_tracethrough = True
stg2.name = 'Stage 2: Cover'
stg2.position = Point(0,0,0)
optics_cover = op_cover_surface(PT, slope_error, specularity_error)
el2 = trapezoidal_secondary_reflector(stg2, optics_cover, receiver_height, receiver_length)

# Stage 3, Heliostats
stg3 = PT.add_stage()
stg3.name = 'Stage 3: Heliostats'
stg3.position = Point(0,0,0)
tilt_list =[]
for p in range(len(panel_positions)):
    optics_heliostat_p = op_heliostat_surface(PT, slope_error, specularity_error, 1, p)

    heliostat_position = [panel_positions[p], 0, panel_height]
    # From the relative positions of the panel to the receiver, the panel's tilt and normal are calculated.
    theta_aim = calculate_theta_aim(Xaim=receiver_position[0], Zaim=receiver_position[2], X0=heliostat_position[0], Z0=heliostat_position[2])
    tilt = calculate_tilt(theta_T, theta_aim)
    tilt_list.append(tilt)
    panel_normal = calculate_panel_normal(tilt)

    el3 = stg3.add_element()
    el3.position = Point(*heliostat_position)
    aim = heliostat_position + 1000*panel_normal
    el3.aim = Point(*aim)
    el3.optic = optics_heliostat_p
    el3.surface_flat()
    el3.aperture_rectangle(panel_width, panel_length)

# Stage 4, Receiver & Secondary Reflector
stg4 = PT.add_stage()
stg4.is_multihit = True
stg4.name = 'Stage 4: Receiver & Secondary Reflector'
stg4.position = Point(0,0,0)
optics_receiver = op_receiver_surface(PT)

el4 = stg4.add_element()
el4.position = Point(*receiver_position)
el4.aim = Point(0,0,0)
el4.optic = optics_receiver
el4.surface_cylindrical(receiver_diameter/2)
el4.aperture_singleax_curve(0, 0, receiver_length) # (inner coordinate of revolved section, outer coordinate of revolved section, 
                                                   # length of revolved section along axis of revolution)
optics_secondary = op_secondaryReflector_surface(PT, slope_error, specularity_error)
el4 = trapezoidal_secondary_reflector(stg4, optics_secondary, receiver_height, receiver_length)

# Simulation Parameters
PT.num_ray_hits = 1e4
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True
PT.is_surface_errors = True
PT.dni= 1000

# When ray data is extracted, the in-built multithreading cannot be used.
PT.run(-1,False)
PT.plot_trace()

# Field Parameters
df = PT.raydata  # Extracting the ray data from the simulation

fictitious_df = df[(df['stage'] == 1) & (df['element'] != 0)]['number'].unique().shape[0]   # Number of rays hitting stage 1
cover_df = df[(df['stage'] == 2) & (df['element'] != 0)]['number'].unique().shape[0]        # Number of rays hitting stage 2
cover_miss = df[(df['stage']==2) & (df['element'] == 0)]['number'].unique().shape[0]        # Number of rays missing stage 2

mirrors_abs = df[(df['stage']==3) & (df['element'] < 0)]['number'].unique().shape[0]     # Number of rays absorbed by stage 3
mirrors_refl = df[(df['stage'] == 3) & (df['element'] > 0)]['number'].unique().shape[0]  # Number of rays reflected by stage 3
mirrors_hits = df[(df['stage']==3) & (df['element'] != 0)]['number'].unique().shape[0]   # Number of rays hitting stage 3
rays_gaps = cover_miss - mirrors_hits # Rays in the space between mirrors

receiver_abs = df[(df['stage'] == 4) & (df['element'] == -1)].shape[0]  # Number of rays hitting receiver from mirrors
receiver_tot = df[(df['stage'] == 4) & (df['element'] != 0)].shape[0] # Number of rays 

# Optical efficiency and power per ray
A_aperture = len(panel_positions) * panel_width * panel_length # Mirrored surface

# Efficiency
ppr_corrected = stg1_width * panel_length * np.cos(theta_T) * PT.dni / (PT.num_ray_hits - rays_gaps)
    
if (mirrors_hits) != 0:
    eta_opt_corrected = receiver_abs * ppr_corrected / (PT.dni * A_aperture)
else:
    eta_opt_corrected = 0   
   
eta_opt_zero = 0.686
eta_opt_rays = receiver_abs / mirrors_hits
IAM = eta_opt_corrected/eta_opt_zero

# Results
print()
print('The sun is at an azimuth of', round(azimuth_deg), 'and a zenith of', round(zenith_deg))
print("Number of rays traced: {:d}".format(PT.raydata.index.size))
print()
print('Number of rays hitting,')
print('Fictitious surface:', fictitious_df)
print('Cover:', cover_df)
print('Mirrors:', mirrors_hits)
print('Receiver:', receiver_abs)
print()
print(f'Rays ratio optical efficiency: {eta_opt_rays*100:.2f}%')
print(f'Corrected power per ray: {ppr_corrected:.2f} W/m2')
print(f'Corrected optical efficiency: {eta_opt_corrected*100:.2f}%')
print(f'Corrected IAM: {IAM*100:.2f}%')