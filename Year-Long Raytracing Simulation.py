
"""
This script implements the soiled results from the csv file soiled_data.csv generated from "Soiled LFR Plant Simulation.py" into
the raytracing simulation set up in LFR Plant Example.py. This allows the efficiency of the plant at each hour of the year to be 
calculated. To speed up this script, not only were many calculations performed in "Soiled LFR Plant Simulation.py" (like finding
the angles of the sun) and then passed to this script in the csv file, but the main loop is also within a multiprocess structure, 
so that multiple cores from the computer's CPU can speed up the calculation. This has resulted in a time reduction of 79%.
"""
# Module Initialisation and Importing CSVs
import numpy as np
import pandas as pd
import time
import os
from multiprocessing import Pool, cpu_count
from pysoltrace import PySolTrace, Point

from raytracing_soiling_functions import calculate_panel_normal
from raytracing_soiling_functions import import_simulation_parameters

from optical_geometrical_setup import op_fictitious_surface
from optical_geometrical_setup import op_cover_surface
from optical_geometrical_setup import op_heliostat_surface
from optical_geometrical_setup import op_receiver_surface
from optical_geometrical_setup import op_secondaryReflector_surface
from optical_geometrical_setup import trapezoidal_secondary_reflector

current_script_path = os.path.abspath(__file__)           
current_directory = os.path.dirname(current_script_path)  
csv_path = os.path.join(current_directory, "CSV Files/soiled_data.csv")
df = pd.read_csv(csv_path)

# LFR Plant Setup 
csv_path = os.path.join(current_directory, "CSV Files/Simulation Parameters.csv")
lat, lon, hour_offset, receiver_height, receiver_length, receiver_diameter, panel_length, panel_width, panel_height, panel_spacing, panels_min_max, slope_error, specularity_error = import_simulation_parameters(pd.read_csv(csv_path))
receiver_position = [0, 0, receiver_height]
panel_positions = np.arange(panels_min_max[0], panels_min_max[1], panel_width + panel_spacing) 
num_heliostats = len(panel_positions)                                

A_aperture = len(panel_positions) * panel_width * panel_length

stg1_length = panel_length 
stg1_width = panel_width * len(panel_positions) + panel_spacing * (len(panel_positions) - 1) 
distance_multiplier = 10                               # Scaling factor which pushes the fictitious surface away from the solar field.
x_shift = (panel_positions[0] + panel_positions[-1])/2 # As the solar field is not exactly centered along the x-axis, there is a shift 
                                                       # required for the aiming algorithm.
z_shift = panel_height/10 # Similarly to the x_shift, there is also a positional shift vertically from the height of the panels.

# Extracting the data required from the soiled_data.csv
azimuths_deg = df['Azimuth [deg]'].to_numpy()
azimuths_rad = np.deg2rad(azimuths_deg)
elevations_deg = df['Elevation [deg]'].to_numpy()
elevations_rad = np.deg2rad(elevations_deg)
num_timesteps = len(df['Date'].to_numpy())
transversal_angles = df['Theta T [rad]'].to_numpy()

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

# Simulating the solar field for every hour of the year
start_time = time.time()
optical_efficiency = []
peak_efficiency = []
efficiency = []

# Main function which is called from the multiprocessing 'Pool', called below. 
# This function is practically the same as the one from LFR Plant Example, where under a specific position of the sun, 
# a raytracing simulation is performed to find the efficiency of the plant.
def ray_trace(i):

    if elevations_deg[i] > 0: # Only simulate raytracing if the sun is above the horizon.
        # Create API class instance
        PT = PySolTrace()

        azimuth_rad = azimuths_rad[i]
        zenith_rad = np.pi/2 - elevations_rad[i]
        theta_T = transversal_angles[i]

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
        el1.position = Point(distance_multiplier*(sun_position[0] + x_shift), 
                             distance_multiplier*sun_position[1], 
                             distance_multiplier*(sun_position[2] + z_shift))
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
        el4 = trapezoidal_secondary_reflector(stg4, optics_secondary, receiver_height, receiver_length)

        # Simulation Parameters
        PT.num_ray_hits = 1e4
        PT.max_rays_traced = PT.num_ray_hits*100
        PT.is_sunshape = True
        PT.is_surface_errors = True
        PT.dni= 1000

        PT.run(-1,False)

        # Field Parameters
        df = PT.raydata  # Extracting the ray data from the simulation

        mirrors_hits = df[(df['stage']==3) & (df['element'] != 0)]['number'].unique().shape[0]   # Number of rays hitting stage 3
        receiver_abs = df[(df['stage'] == 4) & (df['element'] == -1)].shape[0]  # Number of rays hitting receiver from mirrors
        cover_miss = df[(df['stage']==2) & (df['element'] == 0)]['number'].unique().shape[0]        # Number of rays missing stage 2
        rays_gaps = cover_miss - mirrors_hits # Rays in the space between mirrors

        ppr_corrected = stg1_width * panel_length * np.cos(theta_T) * PT.dni / (PT.num_ray_hits - rays_gaps)
    
        if mirrors_hits != 0:
            eta_opt_rays = receiver_abs / mirrors_hits
            eta_opt_corrected = receiver_abs * ppr_corrected / (PT.dni * A_aperture)
        else:
            eta_opt_rays = 0
            eta_opt_corrected = 0

        return [eta_opt_rays, eta_opt_corrected]
    
    else: # When the sun is below the horizon, the efficiency of the plant is 0.
         return [0,0]

csv_data = {}

# Multiprocessing
if __name__ == "__main__":
    cores = cpu_count()  # Number of CPU cores used/ parallel processes. Can be changed to match the appropriate hardware.
    ran = num_timesteps  # Number of inputs the function will process.

    start_multi = time.time()
    with Pool(cores) as pool:
        data_multi = pool.map(ray_trace, range(ran))
    end_multi = time.time()

    print()
    print('Time Taken:', end_multi - start_multi)

    # Appending the results to a csv
    filepath = current_directory + '/CSV Files/raytrace_results.csv'

    uncorrected_efficiency = []
    corrected_efficiency = []
    for i in range(len(data_multi)):
        uncorrected_efficiency.append(data_multi[i][0])
        corrected_efficiency.append(data_multi[i][1])

    csv_data["Uncorrected efficiency"] = uncorrected_efficiency
    csv_data["Corrected efficiency"] = corrected_efficiency

    df = pd.DataFrame(csv_data)
    df.to_csv(filepath, index = False)