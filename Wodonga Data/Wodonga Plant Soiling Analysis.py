#%%
# Initialisation...
"""
Integrating the soiling model from HelioSoil to the LFR plant design in Wodonga.
This script generates a csv file containing information about the tilts of the heliostats along with their respective
soiling. This includes soiling increases for every 5 minutes and the cumulative soiling of the heliostats.
The script's structure follows:
1. Importing all required functions, initialising HelioSoil instances and choosing the plant layout.
2. Calculating a value of hrz0 from the experimental data found in the Training Data folder. 
3. Finding the tilt angles for the heliostats across the entire timeframe.
4. Finding the equivalent soiling area increase from the panel tilts calculated, including a chosen cleaning frequency.
5. Plotting the cumulative soiling of the heliostats across the year.
6. Saving all these results in a csv file.
"""
# Getting the current directory's location
import os
import sys
current_script_path = os.path.abspath(__file__)
wodonga_path = os.path.dirname(current_script_path)
heliosoil_path = wodonga_path.replace(r"\Wodonga Data", "\\HelioSoil")
current_directory = wodonga_path.replace(r"\Wodonga Data", "")

sys.path.append(heliosoil_path)
sys.path.append(current_directory)

# Imported Modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pysolar.solar as solar
import datetime as dt
from datetime import datetime, timezone, timedelta
from copy import deepcopy

import soiling_model.base_models as smb
import soiling_model.fitting as smf
import soiling_model.utilities as smu

from raytracing_soiling_functions import calculate_theta_aim
from raytracing_soiling_functions import calculate_tilt
from raytracing_soiling_functions import import_simulation_parameters
from raytracing_soiling_functions import find_normal

# LFR Plant Setup 
csv_path = os.path.join(wodonga_path, "Wodonga Simulation Parameters.csv")
lat, lon, hour_offset, receiver_height, receiver_length, receiver_diameter, panel_length, panel_width, panel_height, panel_spacing, panels_min_max, slope_error, specularity_error, PSR_divisions, PSR_focal_length, PSR_diameter = import_simulation_parameters(pd.read_csv(csv_path))
receiver_position = [0, 0, receiver_height]
panel_positions = np.arange(panels_min_max[0], panels_min_max[1], panel_width + panel_spacing) 
num_heliostats = len(panel_positions)                                

# File paths to data sheets
file_params = wodonga_path + "\\parameters_wodonga_experiments.xlsx"
file_weather = wodonga_path + "\\Wodonga Data Refactored.xlsx"
wodonga_data = pd.read_excel(file_weather, sheet_name = 'Weather')

Time_data = wodonga_data["Time"]
tz_offset = timezone(timedelta(hours = hour_offset))

def add_timezone(datetime_obj):
    return datetime_obj.replace(tzinfo = tz_offset)

Time_data = Time_data.apply(add_timezone)

# Initialise and load simulation data to the physical model.
sim_data = smb.simulation_inputs(file_weather, dust_type="PM10", verbose=False)
physical_model = smb.physical_base()
physical_model.hrz0 = 50
physical_model.latitude = lat
physical_model.longitude = lon
physical_model.import_site_data_and_constants(file_params)
print('Initialisation Done')
print()

# %%
# Computing hrz0 from data
# The code below is largely taken from the HelioSoil Repository. Visit the HelioSoil repository, and the wodonga_analysis.py script
# for more information. 
reflectometer_incidence_angle = 15       # Angle of incidence of reflectometer
reflectometer_acceptance_angle = 12.5e-3 # Half acceptance angle of reflectance measurements
second_surf = True                       # True if using the second-surface model. Otherwise, use first-surface
d = f"{wodonga_path}/Training Data/"
train_experiments = [0]                  # Indices for training experiements from 0 to len(files) - 1
train_mirrors = ["OE_M1_T00"]            # Which mirrors within the experiemnts are used for
k_factor = None                          # None sets equal to 1.0, "import" imports from the file
dust_type = "PM10"                       

files,all_intervals,exp_mirrors,all_mirrors = smu.get_training_data(d,"experiment_")
orientation = [ [s[1] for s in mirrors] for mirrors in exp_mirrors]

# Feb 2022 (first experiment --- remove last three days after rain started)
all_intervals[0][0] = np.datetime64('2022-02-20T16:20:00')
all_intervals[0][1] = np.datetime64('2022-02-23T17:40:00')

# April 2022 (remove nothing, first/last measurement)
all_intervals[1][0] = np.datetime64('2022-04-21T11:00:00')
all_intervals[1][1] = np.datetime64('2022-04-27T08:30:00')

# Feb 2023 (most recent experiment --- remove very dirty days)
all_intervals[2][0] = np.datetime64('2023-02-09T15:00:00')
all_intervals[2][1] = np.datetime64('2023-02-14T09:45:00')

testing_intervals = all_intervals

Nfiles = len(files)
extract = lambda x,ind: [x[ii] for ii in ind]
files_train = extract(files,train_experiments)
training_intervals = extract(all_intervals,train_experiments)
testing_intervals = list(all_intervals)
t = [t for t in train_experiments]

imodel = smf.semi_physical(file_params)
imodel_constant = smf.constant_mean_deposition(file_params)
sim_data_train = smb.simulation_inputs( files_train,
                                        k_factors=k_factor,
                                        dust_type=dust_type
                                        )
reflect_data_train = smb.reflectance_measurements(  files_train,
                                                    sim_data_train.time,
                                                    number_of_measurements=9.0,
                                                    reflectometer_incidence_angle=reflectometer_incidence_angle,
                                                    reflectometer_acceptance_angle=reflectometer_acceptance_angle,
                                                    import_tilts=True,
                                                    column_names_to_import=train_mirrors
                                                    )
sim_data_train,reflect_data_train = smu.trim_experiment_data(   sim_data_train,
                                                                reflect_data_train,
                                                                training_intervals 
                                                            )
                                                            
sim_data_train,reflect_data_train = smu.trim_experiment_data(   sim_data_train,
                                                                reflect_data_train,
                                                                "reflectance_data" 
                                                            )

imodel.helios_angles(sim_data_train, reflect_data_train, second_surface=second_surf)
imodel.helios.compute_extinction_weights(sim_data_train, imodel.loss_model,
                                         verbose=False, options={'grid_size_x':1000})
ext_weights = imodel.helios.extinction_weighting[0].copy()

imodel_constant.helios_angles(sim_data_train,reflect_data_train,second_surface=second_surf)
file_inds = np.arange(len(files_train))
imodel_constant = smu.set_extinction_coefficients(imodel_constant,ext_weights,file_inds)

log_param_hat,log_param_cov = imodel.fit_mle(   sim_data_train,
                                        reflect_data_train,
                                        transform_to_original_scale=False)

s = np.sqrt(np.diag(log_param_cov))
param_ci = log_param_hat + 1.96*s*np.array([[-1],[1]])
lower_ci = imodel.transform_scale(param_ci[0,:])
upper_ci = imodel.transform_scale(param_ci[1,:])
param_hat = imodel.transform_scale(log_param_hat)
hrz0_mle,sigma_dep_mle = param_hat
print(f'hrz0: {hrz0_mle:.2e} [{lower_ci[0]:.2e}, {upper_ci[0]:.2e}]')
print(f'\sigma_dep: {sigma_dep_mle:.2e} [{lower_ci[1]:.2e},{upper_ci[1]:.2e}] [p.p./day]')

physical_model.hrz0 = hrz0_mle


# %%
# Finding the tilt angles for all the heliostats across the whole time frame
print("Finding all heliostat tilts...")
num_timesteps = len(Time_data)
tilt_angles_rad = np.zeros((num_heliostats, num_timesteps))
incidence_angles_rad = np.zeros((num_heliostats, num_timesteps))
elevation_angles_deg = [] # These are used to collect data, later placed into the csv file.
azimuth_angles_deg = []
transversal_angles = []
for i,date in enumerate(Time_data):

    # Finding the position of the sun
    elevation_deg = solar.get_altitude(lat,lon,date) 
    azimuth_deg = solar.get_azimuth(lat,lon,date)
    elevation_angles_deg.append(elevation_deg)
    azimuth_angles_deg.append(azimuth_deg)

    if elevation_deg > 0:                 # Only finding the tilts of heliostats during the day. 
        zenith_deg = 90 - elevation_deg   # Otherwise, the panels are stowed in the upright position to limit dust deposition.
        elevation_rad, azimuth_rad, zenith_rad = (np.deg2rad(x) for x in [elevation_deg, azimuth_deg, zenith_deg])
        
        # The maths used below are well shown and visualised in the Aiming Strategy for Linear Fresnel Reflectors,
        # in the Helpful Documents folder.
        sun_position = np.array([np.sin(azimuth_rad)*np.sin(zenith_rad), np.cos(azimuth_rad)*np.sin(zenith_rad), np.cos(zenith_rad)])
        sn = sun_position[0:3]/np.linalg.norm(sun_position[0:3])
        theta_T = np.arctan(sn[0]/sn[2])
        transversal_angles.append(theta_T)

        # Finding the tilt of the of the heliostat according to the position of the sun
        for p, x_position in enumerate(panel_positions):
            theta_aim = calculate_theta_aim(Xaim=receiver_position[0], Zaim=receiver_position[2], X0=x_position, Z0=panel_height)
            tilt_angles_rad[p][i] = calculate_tilt(theta_T, theta_aim)
            
            distance_vec = find_normal([x_position, 0, panel_height],[0, 0, receiver_height])
            distance_vec = distance_vec[0:3]/np.linalg.norm(distance_vec[0:3])
            incidence_angles_rad[p][i] = 0.5 * np.arccos(distance_vec.dot(sn))

    else: # When the sun is below the horizon, so it is night and the panels are kept upright.
        for p in range(num_heliostats):
            tilt_angles_rad[p][i] = np.pi/2
            incidence_angles_rad[p][i] = 0
        transversal_angles.append(0)
    
print('All tilt calculations done')

#%%
# Inputting the calculated tilt angles to the physical model
sigma_dep = 0.01 
nominal_reflectance = 1.0

# Using HelioSoil to set up and simulate the soiling of the solar field.
physical_model.helios.tilt = {0: np.rad2deg(tilt_angles_rad)}                                    # Inputting the tilt angles calculated from the entire year.
physical_model.helios.compute_extinction_weights(sim_data, loss_model="geometry", verbose=False) # Computing the optical extinction weights.
physical_model.deposition_flux(sim_data, hrz0=physical_model.hrz0, verbose=False)                # Computing the deposition flux, to populate the pdfqN.
physical_model.calculate_delta_soiled_area(sim_data, sigma_dep=sigma_dep, verbose=False)         # Computing the change in soiled area.
delta_soiled_area = physical_model.helios.delta_soiled_area[0]

#%%
# Finding the equivalent reflectance of each panel for each hour of the year
cleaning_frequency = 14*24 # Hours, representing after how many hours of use the panel's soiled area is reset to zero.

def h(phi):
    return 2/(np.cos(phi))

def calc_reflectance(nominal_reflectance, cumulative_soiled_area, incidence_angle_rad):
    return nominal_reflectance * (1 - cumulative_soiled_area * h(incidence_angle_rad))

# With the cleaning frequency, the cumulative soiled area is calculated.
cumulative_soiled_area = np.zeros_like(delta_soiled_area)
reflectance = np.zeros_like(delta_soiled_area)
for i in range(num_heliostats):
    for t in range(1, num_timesteps):
        cumulative_soiled_area[i,t] = cumulative_soiled_area[i, t-1] + delta_soiled_area[i, t]
        # cumulative_soiled_area[i, t] = 0

        # Constant cleaning Frequency
        # if t % cleaning_frequency == 0:
        #     cumulative_soiled_area[i,t] = 0

        # Conditional Cleaning Frequency
        if (calc_reflectance(nominal_reflectance, cumulative_soiled_area[i,t], incidence_angles_rad[i,t])) < 0.65:
            cumulative_soiled_area[i,t] = 0

        reflectance[i,t] = calc_reflectance(nominal_reflectance, cumulative_soiled_area[i,t], incidence_angles_rad[i,t])
        
for i in range(len(reflectance)):
    reflectance[i][0] = 1.0

#%%
# Plotting the soiling of the heliostats over the year
plt.figure(figsize=(12, 6))
time = np.arange(num_timesteps)

for i in range(num_heliostats):
    plt.plot(time, cumulative_soiled_area[i, :], label=f"Heliostat {i+1}", linewidth=1.5)

plt.xlabel("Time")
plt.ylabel("Cumulative Soiled Area (m²/m²)")
plt.title("Soiling of Heliostats Over Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 6))

for i in range(num_heliostats):
    plt.plot(time, reflectance[i, :], label=f"Heliostat {i+1}", linewidth=1.5)

plt.xlabel("Time")
plt.ylabel("Heliostat Reflectance")
plt.title("Change In Heliostat Reflectance Over Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#%%
# Appending data to a CSV

filepath = wodonga_path + "\\Wodonga Soiled Data.csv"

data = {
    "Date" : Time_data,
    "Azimuth [deg]" : azimuth_angles_deg,
    "Elevation [deg]" : elevation_angles_deg,
    "Theta T [rad]" : transversal_angles,
}

# The 3 for loops below append the heliostat tilts, hourly changes in soiled area, and cumulative soiled area for 
# each heliostat into the data dictionary created above.
for i in range(num_heliostats):
    key = f"Heliostat tilt [deg] {i+1}"
    data[key] = np.rad2deg(tilt_angles_rad[i, :])
for i in range(num_heliostats):
    key = f"Delta Soiled Area [m2/m2] {i+1}"
    data[key] = delta_soiled_area[i]
for i in range(num_heliostats):
    key = f"Cumulative Soiled Area [m2/m2] {i+1}"
    data[key] = cumulative_soiled_area[i]
for i in range(num_heliostats):
    key = f"Heliostat Reflectance {i+1}"
    data[key] = reflectance[i]

df = pd.DataFrame(data)
df.to_csv(filepath, index=False)
print('Done')

#%%