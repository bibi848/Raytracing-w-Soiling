"""
This script graphs all the data collected in both "Soiled LFR Plant Simulation.py" and "Year-Long Raytracing Simulation.py",
and includes some data manipulation, such as calculating the peak efficiency per 24 hour period and the average reflectivies of
the heliostats per hour.
"""


#%% 
# Importing Modules...

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Finding the script's directory to then find the data csv files produced in previous scripts.
current_script_path = os.path.abspath(__file__)           
current_directory = os.path.dirname(current_script_path)  
csv_path = os.path.join(current_directory, "CSV Result Files/soiled_data.csv")
df_soiled_data = pd.read_csv(csv_path)
csv_path = os.path.join(current_directory, "CSV Result Files/raytrace_results.csv")
df_raytrace_results = pd.read_csv(csv_path)

# Extracting the data from the csv file created in Soiled LFR Plant Simulation.py
azimuths_deg = df_soiled_data['Azimuth [deg]'].to_numpy()
azimuths_rad = np.deg2rad(azimuths_deg)
elevations_deg = df_soiled_data['Elevation [deg]'].to_numpy()
elevations_rad = np.deg2rad(elevations_deg)
num_timesteps = len(df_soiled_data['Date'].to_numpy())

num_heliostats = 11

tilt_header_list = []
reflectivity_header_list = []
for i in range(num_heliostats):
    tilt_header = f"Heliostat tilt [deg] {i+1}"
    reflectivity_header = f"Heliostat Reflectivity {i+1}"
    tilt_header_list.append(tilt_header)
    reflectivity_header_list.append(reflectivity_header)

tilts_deg = df_soiled_data[tilt_header_list].to_numpy().T
tilts_rad = np.deg2rad(tilts_deg)
reflectivities = df_soiled_data[reflectivity_header_list].to_numpy().T

# Extracting the data from the csv file created in Year-Long Raytracing Simulation.py
efficiencies = df_raytrace_results["Efficiency"].to_numpy()

#%% 
# Manipulating data...

# As the efficiency has been calculated for each hour of the year, it is broken down into blocks of 24,
# and the peak efficiency is found from each block. This represents the highest efficiency out of the 24 hour period.
days = [efficiencies[i:i + 24] for i in range(0, len(efficiencies), 24)]
peak_efficiencies = []
for day in days:
    peak_efficiencies.append(max(day))

peak_efficiency_times = []
for i in range(len(peak_efficiencies)):
    peak_efficiency_times.append(i*24)

# For each hour of the year, the average reflectivity of the n heliostats is found. 
# This can represent the overall soiling effect on the whole field.
avg_reflectivities = []
for i in range(num_timesteps):
    avg = []
    for p in range(num_heliostats):
        avg.append(reflectivities[p][i])
    avg_reflectivities.append(sum(avg)/len(avg))

#%%
# Plotting data...

plt.figure(figsize=(12, 6))
t = np.arange(num_timesteps)

for i in range(num_heliostats):
    plt.plot(t, reflectivities[i, :], label=f"Heliostat {i+1}", linewidth=1.5)

plt.xlabel("Time (hours)")
plt.ylabel("Cumulative Soiled Area (m²/m²)")
plt.title("Soiling of Heliostats Over Time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


fig, ax1 = plt.subplots(figsize=(12, 6))

ax1.plot(t, avg_reflectivities, color="blue", label="Average Reflectivities")
ax1.set_xlabel("Hours")
ax1.set_ylabel("Average Reflectivities", color="blue")
ax1.tick_params(axis="y", labelcolor="blue")

ax2 = ax1.twinx()  
ax2.plot(peak_efficiency_times, peak_efficiencies, color="red", label="Peak Efficiencies")
ax2.set_ylabel("Peak Efficiency (per day)", color="red")
ax2.tick_params(axis="y", labelcolor="red")

ax1.grid(True)
plt.show()


# %%