'''
This script combines the cleaning optimisation and the raytracing. 
This is because the threshold for cleaning should not be the reflectance, but instead a cleanliness ratio, which can only be found using the optical efficiency.
The efficiency being compared is the DNI weighted efficiency.
'''
#%%
# Getting the current directory's location
import os
import sys
current_script_path = os.path.abspath(__file__)
wodonga_path = os.path.dirname(current_script_path)
heliosoil_path = wodonga_path.replace(r"\Wodonga Data", "\\HelioSoil")
current_directory = wodonga_path.replace(r"\Wodonga Data", "")
filepath = wodonga_path + "\\Clean Optimisation Test 2.csv"

sys.path.append(heliosoil_path)
sys.path.append(current_directory)

# Imported Modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
from multiprocessing import Pool, cpu_count
from pysoltrace import PySolTrace, Point

from raytracing_soiling_functions import import_simulation_parameters
from raytracing_soiling_functions import calculate_panel_normal

from optical_geometrical_setup import op_fictitious_surface
from optical_geometrical_setup import op_cover_surface
from optical_geometrical_setup import op_heliostat_surface
from optical_geometrical_setup import op_receiver_surface
from optical_geometrical_setup import op_secondaryReflector_surface
from optical_geometrical_setup import trapezoidal_secondary_reflector

# LFR Plant Setup 
csv_path = os.path.join(wodonga_path, "Wodonga Simulation Parameters.csv")
lat, lon, hour_offset, receiver_height, receiver_length, receiver_diameter, panel_length, panel_width, panel_height, panel_spacing, panel_positions, slope_error, specularity_error, CPC_depth, aperture_angle = import_simulation_parameters(pd.read_csv(csv_path))
receiver_height -= panel_height   
panel_height = 0
receiver_position = [0, 0, receiver_height]
num_heliostats = len(panel_positions)
aperture = len(panel_positions) * panel_width * panel_length
stg1_length = panel_length 
stg1_width = panel_width * len(panel_positions) + panel_spacing * (len(panel_positions) - 1) 
distance_multiplier = 10                               # Scaling factor which pushes the fictitious surface away from the solar field.

# Extracting the data required from the Wodonga Soiled Data.csv
csv_path = current_directory + "\\Wodonga Data\\Wodonga Soiled Data.csv"
df = pd.read_csv(csv_path)

azimuths_deg = df['Azimuth [deg]'].to_numpy()
azimuths_rad = np.deg2rad(azimuths_deg)
elevations_deg = df['Elevation [deg]'].to_numpy()
elevations_rad = np.deg2rad(elevations_deg)
num_timesteps = len(df['Date'].to_numpy())
transversal_angles = df['Theta T [rad]'].to_numpy()
DNI_values = df['DNI [W/m2]'].to_numpy()
DNI_day_sums = df["DNI Sum Per Day"].to_numpy()

tilt_header_list = []
delta_soiled_area_header_list = []
incident_angle_header_list = []
for i in range(num_heliostats):
    tilt_header = f"Heliostat tilt [deg] {i+1}"
    delta_soiled_area_header = f"Delta Soiled Area [m2/m2] {i+1}"
    incident_angle_header = f"Incident angle [rad] {i+1}"

    tilt_header_list.append(tilt_header)
    delta_soiled_area_header_list.append(delta_soiled_area_header)
    incident_angle_header_list.append(incident_angle_header)

tilts_deg = df[tilt_header_list].to_numpy().T
tilts_rad = np.deg2rad(tilts_deg)
delta_soiled_areas = df[delta_soiled_area_header_list].to_numpy().T
incident_angles_rad = df[incident_angle_header_list].to_numpy().T

csv_path = current_directory + "\\Wodonga Data\\Wodonga Raytrace Results - Clean.csv"
df = pd.read_csv(csv_path)
field_efficiencies_clean = df['Field efficiency'].to_numpy()

# Simulating the solar field for every hour of the year
start_time = time.time()
optical_efficiencies = []
peak_efficiency = []
efficiency = []
cumulative_soiled_area = np.zeros_like(delta_soiled_areas)
cleanliness = np.zeros_like(delta_soiled_areas)
cleanliness_clean = np.ones_like(delta_soiled_areas)

day_soiled = []
day_clean = []

D_soiled_efficiency_full = []
D_clean_efficiency_full = []
field_efficiency_full = []
field_efficiency_clean_full = []

# Cleaning optimisation parameters
cleanliness_ratio = 0.8
nominal_reflectance = 1.0

def shouldClean(soiled_efficiency, clean_efficiency, cleanliness_ratio):
    if (soiled_efficiency/clean_efficiency) < cleanliness_ratio:
        return True
    else:
        return False

def h(phi):
    return 2/(np.cos(phi))

def calc_cleanliness(cumulative_soiled_area, incidence_angle_rad):
    return 1 - cumulative_soiled_area * h(incidence_angle_rad)

def ray_trace(i, cleanlinesses, DNI, sun_position):

    if elevations_deg[i] > 0: # Only simulate raytracing if the sun is above the horizon.
        # Create API class instance
        PT = PySolTrace()
    
        # The XYZ position of the sun is then inputted into the SolTrace API. 
        sun = PT.add_sun()
        sun.position.x = sun_position[0]
        sun.position.y = sun_position[1]
        sun.position.z = sun_position[2]

        # Stage 1, Fictitious Surface
        stg1 = PT.add_stage()
        stg1.is_virtual = True
        stg1.position = Point(0,0,0)

        optics_fictitious = op_fictitious_surface(PT, slope_error, specularity_error)
        el1 = stg1.add_element()
        el1.position = Point(distance_multiplier*sun_position[0], 
                             distance_multiplier*sun_position[1], 
                             distance_multiplier*sun_position[2])
        el1.aim = Point(distance_multiplier*sun_position[0], distance_multiplier*sun_position[1], 0)
        el1.surface_flat()
        el1.aperture_rectangle(stg1_width, stg1_length)
        el1.optic = optics_fictitious

        # Stage 2, Cover
        stg2 = PT.add_stage()
        stg2.is_tracethrough = True
        stg2.position = Point(0,0,0)
        optics_cover = op_cover_surface(PT, slope_error, specularity_error)
        el2 = trapezoidal_secondary_reflector(stg2, optics_cover, receiver_height, receiver_length)

        # Stage 3, Heliostats
        stg3 = PT.add_stage()
        stg3.position = Point(0,0,0)
        for p in range(num_heliostats):
            optics_heliostat_p = op_heliostat_surface(PT, slope_error, specularity_error, cleanlinesses[p][i], p)

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
        PT.num_ray_hits = 1e3
        PT.max_rays_traced = PT.num_ray_hits*100
        PT.is_sunshape = True
        PT.is_surface_errors = True
        PT.dni= DNI

        PT.run(-1,False)

        # Field Parameters
        df = PT.raydata  # Extracting the ray data from the simulation

        mirrors_hits = df[(df['stage']==3) & (df['element'] != 0)]['number'].unique().shape[0] # Number of rays hitting stage 3
        receiver_abs = df[(df['stage'] == 4) & (df['element'] == -1)].shape[0]                 # Number of rays hitting receiver from mirrors
    
        if mirrors_hits != 0:
            ppr = PT.powerperray
            field_efficiency = (receiver_abs * ppr) / (PT.dni * aperture)

        else:
            field_efficiency = 0

        return field_efficiency
    
    else: # When the sun is below the horizon, the efficiency of the plant is 0.
         return 0
    
def multi_func(i):

    global day_clean
    global day_soiled

    # Reflectance calculation
    for p in range(num_heliostats):
        if i != 0:
            cumulative_soiled_area[p, i] = cumulative_soiled_area[p, i-1] + delta_soiled_areas[p, i]
        cleanliness[p, i] = calc_cleanliness(cumulative_soiled_area[p, i], incident_angles_rad[p, i])

    # Raytracing to find optical efficiency
    azimuth_rad = azimuths_rad[i]
    zenith_rad = np.pi/2 - elevations_rad[i]
    DNI = DNI_values[i]
    sun_position = np.array([np.sin(azimuth_rad)*np.sin(zenith_rad), np.cos(azimuth_rad)*np.sin(zenith_rad), np.cos(zenith_rad)])

    field_efficiency = ray_trace(i, cleanliness, DNI, sun_position)
    field_efficiency_clean = field_efficiencies_clean[i]

    D_optical_efficiency = field_efficiency * DNI
    D_optical_efficiency_clean = field_efficiency_clean * DNI

    field_efficiency_full.append(field_efficiency)
    field_efficiency_clean_full.append(field_efficiency_clean)
    
    day_soiled.append(D_optical_efficiency)
    day_clean.append(D_optical_efficiency_clean)

    # if DNI != 0:
    #     soiled_efficiency = D_optical_efficiency / DNI
    #     clean_efficiency = D_optical_efficiency_clean / DNI

    # if shouldClean(soiled_efficiency, clean_efficiency, cleanliness_ratio):
    #         cumulative_soiled_area[:, i] = 0
    
    # return [soiled_efficiency, clean_efficiency, optical_efficiency, optical_efficiency_clean]

    # Check cleanliness
    if i % 288 == 0 and i != 0:
        
        soiled_efficiency = sum(day_soiled) / DNI_day_sums[int(i/288)-1]
        clean_efficiency = sum(day_clean) / DNI_day_sums[int(i/288)-1]

        if shouldClean(soiled_efficiency, clean_efficiency, cleanliness_ratio):
            cumulative_soiled_area[:, i] = 0

        day_clean = []
        day_soiled = []

        D_soiled_efficiency_full.append(soiled_efficiency)
        D_clean_efficiency_full.append(clean_efficiency)


for i in range(num_timesteps):
    print(i)
    multi_func(i)

print()
print('Inputting into csv')

extend = len(field_efficiency_full) - len(D_soiled_efficiency_full)
D_soiled_efficiency_full.extend([0] * extend)
D_clean_efficiency_full.extend([0] * extend)

data = {
    "Optical Efficiency" : field_efficiency_full,
    "Optical Efficiency - Clean" : field_efficiency_clean_full,
    "DNI Corrected Efficiency [per day]" : D_soiled_efficiency_full,
    "DNI Corrected Efficiency (clean) [per day]" : D_clean_efficiency_full
}

for i in range(num_heliostats):
    key = f"Cumulative Soiled Area [m2/m2] {i+1}"
    data[key] = cumulative_soiled_area[i]
for i in range(num_heliostats):
    key = f"Heliostat Cleanliness {i+1}"
    data[key] = cleanliness[i]


df = pd.DataFrame(data)
df.to_csv(filepath, index = False)


# data = {}

# # Multiprocessing
# if __name__ == "__main__":
#     cores = cpu_count()  # Number of CPU cores used/ parallel processes. Can be changed to match the appropriate hardware.
#     ran = num_timesteps  # Number of inputs the function will process.

#     start_multi = time.time()
#     with Pool(cores) as pool:
#         data_multi = pool.map(multi_func, range(ran))
#     end_multi = time.time()

#     print()
#     print('Time Taken:', end_multi - start_multi)

#     DNI_efficiency = []
#     DNI_efficiency_clean = []
#     optical_efficiency = []
#     optical_efficiency_clean = []
#     for i in range(len(data_multi)):
#         DNI_efficiency.append(data_multi[i][0])
#         DNI_efficiency_clean.append(data_multi[i][1])
#         optical_efficiency.append(data_multi[i][2])
#         optical_efficiency_clean.append(data_multi[i][3])

#     data = {
#         "Optical Efficiency" : optical_efficiency,
#         "Optical Efficiency (clean)" : optical_efficiency_clean,
#         "DNI Corrected Efficiency" : DNI_efficiency,
#         "DNI Corrected Efficiency (clean)" : DNI_efficiency_clean
#     }

#     for i in range(num_heliostats):
#         key = f"Heliostat Reflectance {i+1}"
#         data[key] = reflectance[i]

#     df = pd.DataFrame(data)
#     df.to_csv(filepath, index = False)
