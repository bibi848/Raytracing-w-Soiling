#%%
"""
Integrating the soiling model from HelioSoil to the LFR plant design.
This script attempts to model the efficiency of the plant across the whole year (2018) from the data collected in HelioSoil.
The script's structure follows:
1. Importing all required functions and initialising the models used in HelioSoil.
2. Calculating all the Fresnels' tilts across the year, generating data which represents the soiling across the whole year.
3. Finding the reflectivity of the panels for each hour of the year.
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
import pysolar.solar as solar
from pysoltrace import PySolTrace, Point
import datetime as dt
from datetime import datetime, timedelta
import pysolar.radiation as radiation
import warnings
warnings.filterwarnings("ignore", message="no explicit representation of timezones available for np.datetime64")

import soiling_model.base_models as smb

from raytracing_soiling_functions import calculate_theta_aim
from raytracing_soiling_functions import calculate_tilt

# LFR Plant Setup
lat, lon = -31.2,136.816667

receiver_height = 4.5    # [m]
receiver_length = 12     # [m]
receiver_diameter = 0.15 # [m]
receiver_position = [0, 0, receiver_height]

panel_length = 12        # [m]
panel_width = 0.5        # [m]
panel_height = 0.5       # [m]
panel_spacing = 0.2      # [m]
panel_positions = np.arange(-3.5, 3.75, panel_width + panel_spacing) # Describes the x coordinate for each of the mirros, 
num_heliostats = len(panel_positions)                                # ranging from -3.5m to 3.75m.


# Time parameters
timezoneOffset = dt.timedelta(hours = 9.5)
start_date = datetime(2018, 1, 1, hour=0, minute=0, second=0, tzinfo=dt.timezone(timezoneOffset))
end_date = datetime(2018, 12, 31, hour=23, minute=0, second=0, tzinfo=dt.timezone(timezoneOffset))
datetime_list = []
date = start_date
while date <= end_date:
    datetime_list.append(date)
    date += timedelta(hours=1)
num_timesteps = len(datetime_list)
print('Datetime list done')

# File paths to data sheets
file_params = heliosoil_path + "/woomera_demo/parameters.xlsx"
file_weather = heliosoil_path + "/woomera_demo/woomera_data.xlsx"

# Initialise and load simulation data to the physical model.
sim_data = smb.simulation_inputs(file_weather, dust_type="PM10", verbose=False)
physical_model = smb.physical_base()
physical_model.import_site_data_and_constants(file_params)

# Finding the tilt angles for all heliostats across the whole year
print('Finding all heliostat tilts...')
tilt_angles_rad = np.zeros((num_heliostats, num_timesteps))
for i,date in enumerate(datetime_list):

    # Finding the position of the sun
    elevation_deg = solar.get_altitude(lat,lon,date) 

    if elevation_deg > 0:
        azimuth_deg = solar.get_azimuth(lat,lon,date) 
        zenith_deg = 90 - elevation_deg  
        elevation_rad, azimuth_rad, zenith_rad = (np.deg2rad(x) for x in [elevation_deg, azimuth_deg, zenith_deg])

        sun_position = np.array([np.sin(azimuth_rad)*np.sin(zenith_rad), np.cos(azimuth_rad)*np.sin(zenith_rad), np.cos(zenith_rad)])
        sn = sun_position[0:3]/np.linalg.norm(sun_position[0:3])
        theta_T = np.arctan(sn[0]/sn[2])

        # Finding the tilt of the of the heliostat according to the position of the sun
        for p, x_position in enumerate(panel_positions):
            theta_aim = calculate_theta_aim(Xaim=receiver_position[0], Zaim=receiver_position[2], X0=x_position, Z0=panel_height)
            tilt_angles_rad[p][i] = calculate_tilt(theta_T, theta_aim)
    else:
        for p, x_position in enumerate(panel_positions):
            tilt_angles_rad[p][i] = np.pi/2
    
print('All tilt calculations done')

#%%
# Finding the equivalent reflectivity of each panel for each hour of the year
panel_reflectivities = np.zeros((num_heliostats, num_timesteps))
sigma_dep = 0.01

physical_model.helios.tilt = {0: np.rad2deg(tilt_angles_rad)}
physical_model.helios.compute_extinction_weights(sim_data, loss_model="geometry", verbose=True)
physical_model.deposition_flux(sim_data, hrz0=physical_model.hrz0, verbose=True)
physical_model.calculate_delta_soiled_area(sim_data, sigma_dep=sigma_dep, verbose=True)
delta_soiled_area = physical_model.helios.delta_soiled_area[0]
cumulative_soiled_area = np.cumsum(delta_soiled_area, axis=1)

#%%
# Plotting the soiling of the heliostats over the year
time = np.arange(num_timesteps)

# Plot the delta_soiled_area for each heliostat
plt.figure(figsize=(12, 8))

for i in range(num_heliostats):
    plt.plot(time, cumulative_soiled_area[i, :], label=f"Heliostat {i+1}", linewidth=1.5)

# Add labels and title
plt.xlabel("Time (hours)")
plt.ylabel("Soiled Area (m²/m²)")
plt.title("Soiling of Heliostats Over Time")
plt.legend(loc='upper left')
plt.grid(True)

# Show the plot
plt.tight_layout()
plt.show()

#%%

