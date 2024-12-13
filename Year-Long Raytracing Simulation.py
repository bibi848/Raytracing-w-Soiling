#%%
# Module Initialisation and Importing CSV
"""
This script uses the CSV file generated in "Soiled LFR Plant Simulation.py" to take the tilt angles and the reflectivity values
and applies it to the raytracing calculations. 
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time as tmr
from pysoltrace import PySolTrace, Point
import warnings
warnings.filterwarnings("ignore", message="no explicit representation of timezones available for np.datetime64")

from raytracing_soiling_functions import calculate_theta_aim
from raytracing_soiling_functions import calculate_tilt
from raytracing_soiling_functions import calculate_panel_normal

from optical_geometrical_setup import op_fictitious_surface
from optical_geometrical_setup import op_cover_surface
from optical_geometrical_setup import op_heliostat_surface
from optical_geometrical_setup import op_receiver_surface
from optical_geometrical_setup import op_secondaryReflector_surface
from optical_geometrical_setup import trapezoidal_secondary_reflector

import os
current_script_path = os.path.abspath(__file__)           
current_directory = os.path.dirname(current_script_path)  
csv_path = os.path.join(current_directory, "CSV Result Files/soiled_data.csv")
df = pd.read_csv(csv_path)

# Plant layout
lat, lon = -31.2,136.816667 # Woomera
receiver_height = 4.5    # [m]
receiver_length = 12     # [m]
receiver_diameter = 0.15 # [m]
receiver_position = [0, 0, receiver_height]

panel_length = 12        # [m]
panel_width = 0.5        # [m]
panel_height = 0         # [m], in reality this is innaccurate. However, when aiming the fictitious surface this is necessary.
panel_spacing = 0.2      # [m]

panel_positions = np.arange(-3.5, 3.75, panel_width + panel_spacing) # Describes the x-coordinate for each of the mirrors, 
                                                                     # ranging from -3.5 [m] to 3.75 [m].
num_heliostats = len(panel_positions)
slope_error = 0.1        # [mrad]
specularity_error = 0.1  # [mrad]

stg1_length = panel_length 
stg1_width = panel_width * len(panel_positions) + panel_spacing * (len(panel_positions) - 1) 
distance_multiplier = 10                               # Scaling factor which pushes the fictitious surface away from the solar field.
x_shift = (panel_positions[0] + panel_positions[-1])/2 # As the solar field is not exactly centered along the x-axis, there is a shift 
                                                       # required for the aiming algorithm.

azimuths_deg = df['Azimuth [deg]'].to_numpy()
azimuths_rad = np.deg2rad(azimuths_deg)
elevations_deg = df['Elevation [deg]'].to_numpy()
elevations_rad = np.deg2rad(elevations_deg)
num_timesteps = len(df['Date'].to_numpy())

tilt_header_list = []
reflectivity_header_list = []
for i in range(num_heliostats):
    tilt_header = f"Heliostat tilt [deg] {i+1}"
    reflectivity_header = f"Heliostat Reflectivity {i+1}"
    tilt_header_list.append(tilt_header)
    reflectivity_header_list.append(reflectivity_header)

tilts_deg = df[tilt_header_list].to_numpy().T
tilts_rad = np.deg2rad(tilts_deg)
reflectivities = df[reflectivity_header_list].to_numpy().T

#%%
# Simulating the solar field for every hour of the year
start_time = tmr.time()
optical_efficiency = []
peak_efficiency = []
efficiency = []

for i in range(num_timesteps):

    print(i)

    if (i % 24 == 0) and (i > 0):
        peak_efficiency.append(max(efficiency))
        efficiency = []

    if elevations_deg[i] > 0:
        # Create API class instance
        PT = PySolTrace()

        azimuth_rad = azimuths_rad[i]
        zenith_rad = np.pi/2 - elevations_rad[i]
    
        # Describing the sun's position in terms of the azimuth and zenith. The full breakdown for this result is shown in
        # Aiming Strategy for LFRs document. 
        sun_position = np.array([np.sin(azimuth_rad)*np.sin(zenith_rad), np.cos(azimuth_rad)*np.sin(zenith_rad), np.cos(zenith_rad)])
    
        # The XYZ position of the sun is then inputted into the SolTrace API. 
        sun = PT.add_sun()
        sun.position.x = sun_position[0]
        sun.position.y = sun_position[1]
        sun.position.z = sun_position[2]

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
        for p in range(num_heliostats):
            optics_heliostat_p = op_heliostat_surface(PT, slope_error, specularity_error, reflectivities[p][i], p)

            heliostat_position = [panel_positions[p], 0, panel_height]
            panel_normal = calculate_panel_normal(tilts_rad[p][i])

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
        el4 = trapezoidal_secondary_reflector(stg3, optics_secondary, receiver_height, receiver_length)

        # Simulation Parameters
        PT.num_ray_hits = 1e4
        PT.max_rays_traced = PT.num_ray_hits*100
        PT.is_sunshape = True
        PT.is_surface_errors = True
        PT.dni= 1000

        # When ray data is extracted, multithreading cannot be used.
        PT.run(-1,False)

        # Field Parameters
        df = PT.raydata  # Extracting the ray data from the simulation

        mirrors_hits = df[(df['stage']==3) & (df['element'] != 0)]['number'].unique().shape[0]   # Number of rays hitting stage 3
        receiver_abs = df[(df['stage'] == 4) & (df['element'] == -1)].shape[0]  # Number of rays hitting receiver from mirrors
    
        if mirrors_hits != 0:
            eta_opt_rays = receiver_abs / mirrors_hits
        else:
            eta_opt_rays = 0

        #optical_efficiency.append(eta_opt_rays)
        efficiency.append(eta_opt_rays)
    
    else:
        #optical_efficiency.append(0)
        efficiency.append(0)

end_time = tmr.time()
print('Simulation Took', (end_time - start_time), 'seconds')

#%%
# Plotting Collected Data

time = np.arange(num_timesteps)
ax2_x_axis = np.arange(len(peak_efficiency))
fig, ax1 = plt.subplots()

for i in range(num_heliostats):
    ax1.plot(time, reflectivities[i, :], label=f"Heliostat {i+1}", linewidth=1.5)
ax1.set_xlabel('Time [Hours]')
ax1.set_ylabel('Heliostat Reflectivities', color = 'blue')
ax1.tick_params(axis='y', labelcolor = 'blue')

ax2 = ax1.twinx()
ax2.plot(ax2_x_axis, peak_efficiency, color = 'red', label = 'Peak Efficiency')
ax2.set_ylabel('Peak Efficiency per 24 hour period', color = 'red')
ax2.tick_params(axis='y', labelcolor='red')

ax1.grid(True)
plt.show()

# %%
