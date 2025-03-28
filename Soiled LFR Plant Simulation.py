#%%
# Initialisation...
"""
Integrating the soiling model from HelioSoil to the LFR plant design.
This script generates a csv file containing information about the tilts of the heliostats across 2018 along with their respective
soiling. This includes hourly soiling increases and the cumulative soiling of the heliostats.
The script's structure follows:
1. Importing all required functions, initialising HelioSoil instances and choosing the plant layout.
2. Generating the datetime list for the year 2018. The timestep used is every hour (8760 hours in total).
3. Finding the tilt angles for the heliostats across the entire year.
4. Finding the equivalent soiling area increase from the panel tilts calculated, including a chosen cleaning frequency.
5. Plotting the cumulative soiling of the heliostats across the year.
6. Saving all these results in a csv file.
"""

# Getting the current directory's location
import os
import sys
current_script_path = os.path.abspath(__file__)           
current_directory = os.path.dirname(current_script_path)  

heliosoil_path = os.path.join(current_directory, "HelioSoil")
sys.path.append(heliosoil_path)

# Imported Modules
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pysolar.solar as solar
import datetime as dt
from datetime import datetime, timedelta

from pvlib.location import Location
from pvlib.clearsky import ineichen

import soiling_model.base_models as smb

from raytracing_soiling_functions import calculate_theta_aim
from raytracing_soiling_functions import calculate_tilt
from raytracing_soiling_functions import import_simulation_parameters
from raytracing_soiling_functions import find_normal

# LFR Plant Setup 
csv_path = os.path.join(current_directory, "CSV Files/Simulation Parameters.csv")
lat, lon, hour_offset, receiver_height, receiver_length, receiver_diameter, panel_length, panel_width, panel_height, panel_spacing, number_of_modules, panels_per_module, slope_error, specularity_error = import_simulation_parameters(pd.read_csv(csv_path))
receiver_height -= panel_height   
panel_height = 0
num_heliostats = number_of_modules * panels_per_module
stg1_width = panel_width * num_heliostats + panel_spacing * (num_heliostats - 1) 

panel_x_shift = panel_width + panel_spacing
panel_positions = []
panel_pos = -stg1_width/2 + panel_width/2

for i in range(num_heliostats):
    panel_positions.append(panel_pos)
    panel_pos += panel_x_shift

panel_positions_b = [panel_positions[i:i + panels_per_module] for i in range(0, len(panel_positions), panels_per_module)]

receiver_positions = []
for module in panel_positions_b:
    avg = (module[0] + module[-1]) / 2
    receiver_positions.append(avg)

# File paths to data sheets
file_params = heliosoil_path + "/woomera_demo/parameters.xlsx"
file_weather = heliosoil_path + "/woomera_demo/woomera_data.xlsx"

# Initialise and load simulation data to the physical model.
sim_data = smb.simulation_inputs(file_weather, dust_type="PM10", verbose=False)
physical_model = smb.physical_base()
physical_model.import_site_data_and_constants(file_params)
print('Initialisation Done')
print()

#%%
# Generating the datetime list
# Adding datetime objects every hour between the start and end date to datetime_list
print('Generating Datetime list...')
timezoneOffset = dt.timedelta(hours = hour_offset)
start_date = datetime(2018, 1, 1, hour=0, minute=0, second=0, tzinfo=dt.timezone(timezoneOffset))
end_date = datetime(2018, 12, 31, hour=23, minute=0, second=0, tzinfo=dt.timezone(timezoneOffset))
datetime_list = []
date = start_date
while date <= end_date:
    datetime_list.append(date)
    date += timedelta(hours=1)
num_timesteps = len(datetime_list)
print('Datetime list done')

#%%
# Finding the DNI for every hour of the year
location = Location(lat, lon, dt.timezone(timezoneOffset), 0)
times = pd.date_range(start='2018-01-01', end='2019-01-01', freq='1h', tz=dt.timezone(timezoneOffset))
clear_sky = location.get_clearsky(times, model = 'ineichen')
DNI = list(clear_sky['dni'])
del DNI[-1]

#%%
# Finding the tilt angles for all heliostats across the whole year
print('Finding all heliostat tilts...')
tilt_angles_rad = np.zeros((num_heliostats, num_timesteps))
incidence_angles_rad = np.zeros((num_heliostats, num_timesteps))
elevation_angles_deg = [] # These are used to collect data, later placed into the csv file.
azimuth_angles_deg = []
transversal_angles = []
for i,date in enumerate(datetime_list):

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
        q = 0
        for p in range(num_heliostats):
            if p != 0 and p != 1 and panels_per_module % p == 0:
                q += 1
        
            receiver_position = [receiver_positions[q], 0, receiver_height]
            x_position = panel_positions[p]

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
# Finding the equivalent reflectance of each panel for each hour of the year
sigma_dep = 0.01 
nominal_reflectance = 1.0

# Using HelioSoil to set up and simulate the soiling of the solar field.
physical_model.helios.tilt = {0: np.rad2deg(tilt_angles_rad)}                                    # Inputting the tilt angles calculated from the entire year.
physical_model.helios.compute_extinction_weights(sim_data, loss_model="geometry", verbose=False) # Computing the optical extinction weights.
physical_model.deposition_flux(sim_data, hrz0=physical_model.hrz0, verbose=False)                # Computing the deposition flux, to populate the pdfqN.
physical_model.calculate_delta_soiled_area(sim_data, sigma_dep=sigma_dep, verbose=False)         # Computing the change in soiled area.
delta_soiled_area = physical_model.helios.delta_soiled_area[0]

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

        # # Constant cleaning Frequency
        # if t % cleaning_frequency == 0:
        #     cumulative_soiled_area[i,t] = 0

        # Conditional Cleaning Frequency
        if (calc_reflectance(nominal_reflectance, cumulative_soiled_area[i,t], incidence_angles_rad[i,t])) < 0.7:
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

plt.xlabel("Time (hours)")
plt.ylabel("Cumulative Soiled Area (m²/m²)")
plt.title("Soiling of Heliostats Over Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 6))

for i in range(num_heliostats):
    plt.plot(time, reflectance[i, :], label=f"Heliostat {i+1}", linewidth=1.5)

plt.xlabel("Time (hours)")
plt.ylabel("Heliostat Reflectance")
plt.title("Change In Heliostat Reflectance Over Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

#%%
# Appending data to a CSV

filepath = current_directory + '/CSV Files/soiled_data.csv'

data = {
    "Date" : datetime_list,
    "Azimuth [deg]" : azimuth_angles_deg,
    "Elevation [deg]" : elevation_angles_deg,
    "DNI [W/m2]" : DNI,
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

# %%
