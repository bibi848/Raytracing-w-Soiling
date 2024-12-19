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
csv_path = os.path.join(current_directory, "CSV Files/soiled_data.csv")
df_soiled_data = pd.read_csv(csv_path)
csv_path = os.path.join(current_directory, "CSV Files/raytrace_results.csv")
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
uncorrected_efficiencies = df_raytrace_results["Uncorrected efficiency"].to_numpy()
corrected_efficiencies = df_raytrace_results["Corrected efficiency"].to_numpy()

#%% 
# Manipulating data...

# As the efficiency has been calculated for each hour of the year, it is broken down into blocks of 24,
# and the peak efficiency is found from each block. This represents the highest efficiency out of the 24 hour period.
days = [uncorrected_efficiencies[i:i + 24] for i in range(0, len(uncorrected_efficiencies), 24)]
peak_efficiencies = []
avg_efficiencies = []
for day in days:
    peak_efficiencies.append(max(day))
    avg_efficiencies.append(np.mean(day))

peak_efficiency_times = []
for i in range(len(peak_efficiencies)):
    peak_efficiency_times.append(i*24)

days_corrected = [corrected_efficiencies[i:i + 24] for i in range(0, len(corrected_efficiencies), 24)]
peak_c_efficiencies = []
avg_c_efficiencies = []
for day in days_corrected:
    peak_c_efficiencies.append(max(day))
    avg_c_efficiencies.append(np.mean(day))

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

# Plotting the change in heliostat reflectivity over the year
plt.figure(figsize=(12, 6))
t = np.arange(num_timesteps)

for i in range(num_heliostats):
    plt.plot(t, reflectivities[i, :], label=f"Heliostat {i+1}", linewidth=1.5)

plt.xlabel("Time (hours)")
plt.ylabel("Reflectivity")
plt.title("Change in Heliostat Reflectivities")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Overlaying the change in reflectivity with the uncorrected peak efficiencies

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t, avg_reflectivities, color="blue", label="Average Reflectivities")
ax1.set_xlabel("Hours")
ax1.set_ylabel("Average Reflectivity", color="blue")
ax1.set_title('Comparing the Average Field Reflectivity to the Peak Efficiencies')
ax1.tick_params(axis="y", labelcolor="blue")

ax2 = ax1.twinx()  
ax2.plot(peak_efficiency_times, peak_efficiencies, color="red", label="Uncorrected")
ax2.plot(peak_efficiency_times, peak_c_efficiencies, color="green", label="Corrected")
ax2.set_ylabel("Peak Efficiency (per day)", color="black")
ax2.tick_params(axis="y", labelcolor="black")
ax2.legend(loc='lower right')
ax1.grid(True)
plt.show()

# Plant Efficiency Over a 24 Hour Period
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t[:24], uncorrected_efficiencies[:24], color="red", label = 'Uncorrected')
ax1.plot(t[:24], corrected_efficiencies[:24], color="green", label = 'Corrected')
ax1.set_xlabel("Hours")
ax1.set_ylabel("Efficiency Across 24 Hour Period", color="black")
ax1.set_title('Plant Efficiency Over a 24 Hour Period')
ax1.tick_params(axis="y", labelcolor="black")
ax1.legend()
ax1.grid(True)
plt.show()

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t[2812:2851], uncorrected_efficiencies[2812:2851], color="red", label = 'Uncorrected')
ax1.plot(t[2812:2851], corrected_efficiencies[2812:2851], color="green", label = 'Corrected')
ax1.set_xlabel("Hours")
ax1.set_ylabel("Test", color="black")
ax1.set_title('test')
ax1.tick_params(axis="y", labelcolor="black")
ax1.legend()
ax1.grid(True)
plt.show()

# Average efficiencies instead of peak
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t, avg_reflectivities, color="blue", label="Average Reflectivities")
ax1.set_xlabel("Hours")
ax1.set_ylabel("Average Reflectivity", color="blue")
ax1.set_title('Comparing the Average Field Reflectivity to the Avg Efficiencies')
ax1.tick_params(axis="y", labelcolor="blue")

ax2 = ax1.twinx()  
ax2.plot(peak_efficiency_times, avg_efficiencies, color="red", label="Uncorrected")
ax2.plot(peak_efficiency_times, avg_c_efficiencies, color="green", label="Corrected")
ax2.set_ylabel("Average Efficiency (per 24 hour period)", color="black")
ax2.tick_params(axis="y", labelcolor="black")
ax2.legend(loc='lower right')
ax1.grid(True)
plt.show()
# %%
