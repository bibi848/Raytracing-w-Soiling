"""
This script takes the data produced both by the ray tracing and soiling analysis scripts and provides a space where the data can
be analysed and visulised. 
"""
#%% 
# Importing Modules...

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Finding the script's directory to then find the data csv files produced in previous scripts.
current_script_path = os.path.abspath(__file__)           
current_directory = os.path.dirname(current_script_path)[:50] 
csv_path = os.path.join(current_directory, "Wodonga Data\\Wodonga Soiled Data.csv")
df_soiled_data = pd.read_csv(csv_path)
csv_path = os.path.join(current_directory, "Wodonga Data\\Wodonga Clean Data.csv")
df_clean_data = pd.read_csv(csv_path)
csv_path = os.path.join(current_directory, "Wodonga Data\\Wodonga Raytrace Results.csv")
df_raytrace_results = pd.read_csv(csv_path)
csv_path = os.path.join(current_directory, "Wodonga Data\\Wodonga Raytrace Results - Clean.csv")
df_raytrace_results_clean = pd.read_csv(csv_path)

# Extracting the data from the csv file created in Wodonga Plant Soiling Analysis.py
azimuths_deg = df_soiled_data['Azimuth [deg]'].to_numpy()
azimuths_rad = np.deg2rad(azimuths_deg)
elevations_deg = df_soiled_data['Elevation [deg]'].to_numpy()
elevations_rad = np.deg2rad(elevations_deg)
num_timesteps = len(df_soiled_data['Date'].to_numpy())

num_heliostats = 11

tilt_header_list = []
reflectance_header_list = []
for i in range(num_heliostats):
    tilt_header = f"Heliostat tilt [deg] {i+1}"
    reflectance_header = f"Heliostat Reflectance {i+1}"
    tilt_header_list.append(tilt_header)
    reflectance_header_list.append(reflectance_header)

tilts_deg = df_soiled_data[tilt_header_list].to_numpy().T
reflectances = df_soiled_data[reflectance_header_list].to_numpy().T

# Extracting the data from the csv file created in Wodonga Ray Tracing.py
uncorrected_efficiencies = df_raytrace_results["Uncorrected efficiency"].to_numpy()
corrected_efficiencies = df_raytrace_results["Corrected efficiency"].to_numpy()
uncorrected_efficiencies_clean = df_raytrace_results_clean["Uncorrected efficiency"].to_numpy()
corrected_efficiencies_clean = df_raytrace_results_clean["Corrected efficiency"].to_numpy()

# Removing all zero efficiency values (night-time)
uncorrected_efficiencies_no_zeroes = uncorrected_efficiencies[uncorrected_efficiencies != 0]
corrected_efficiencies_no_zeroes = corrected_efficiencies[corrected_efficiencies != 0]
corrected_efficiencies_clean_no_zeroes = corrected_efficiencies_clean[corrected_efficiencies_clean != 0]
uncorrected_efficiencies_clean_no_zeroes = uncorrected_efficiencies_clean[uncorrected_efficiencies_clean != 0]
#%% 
# Manipulating data...

# Average Field Reflectance
avg_reflectances = []
for i in range(num_timesteps):
    avg = []
    for p in range(num_heliostats):
        avg.append(reflectances[p][i])
    avg_reflectances.append(sum(avg) / len(avg))

# As the efficiency has been calculated for every 5 minutes, it is broken down into blocks of 288, which represent 24 hour periods.
efficiency_data_length = len(uncorrected_efficiencies)
efficiency_data_length_no_zeroes = len(uncorrected_efficiencies_no_zeroes)
full_day_sample_number = (24 * 60) / 5  # Number of minutes in a full 24 hour period, divided by the sample time step (5 minutes)
x_shift = num_timesteps / int(full_day_sample_number)

c_efficiencies_full = []
uc_efficiencies_full = []
c_efficiencies_full_clean = []
uc_efficiencies_full_clean = []

c_efficiencies_day = []
uc_efficiencies_day = []
c_efficiencies_day_clean = []
uc_efficiencies_day_clean = []

for i in range(len(df_raytrace_results["Corrected efficiency"])):

    c_efficiency = df_raytrace_results["Corrected efficiency"][i]
    uc_efficiency = df_raytrace_results["Uncorrected efficiency"][i]
    c_efficiency_clean = df_raytrace_results_clean["Corrected efficiency"][i]
    uc_efficiency_clean = df_raytrace_results_clean["Uncorrected efficiency"][i]

    if c_efficiency != 0:

        c_efficiencies_day.append(c_efficiency)
        uc_efficiencies_day.append(uc_efficiency)
        c_efficiencies_day_clean.append(c_efficiency_clean)
        uc_efficiencies_day_clean.append(uc_efficiency_clean)

        if i < len(df_raytrace_results["Corrected efficiency"])-1 and df_raytrace_results["Corrected efficiency"][i+1] < 0.01:

            c_efficiencies_full.append(c_efficiencies_day)
            uc_efficiencies_full.append(uc_efficiencies_day)
            c_efficiencies_full_clean.append(c_efficiencies_day_clean)
            uc_efficiencies_full_clean.append(uc_efficiencies_day_clean)

            c_efficiencies_day = []
            uc_efficiencies_day = []
            c_efficiencies_day_clean = []
            uc_efficiencies_day_clean = []

c_efficiencies_avg = []
uc_efficiencies_avg = []
c_efficiencies_avg_clean = []
uc_efficiencies_avg_clean = []

c_efficiencies_peak = []
uc_efficiencies_peak = []
c_efficiencies_peak_clean = []
uc_efficiencies_peak_clean = []

for i in range(len(c_efficiencies_full)):
    c_efficiencies_avg.append(np.mean(c_efficiencies_full[i]))
    uc_efficiencies_avg.append(np.mean(uc_efficiencies_full[i]))
    c_efficiencies_avg_clean.append(np.mean(c_efficiencies_full_clean[i]))
    uc_efficiencies_avg_clean.append(np.mean(uc_efficiencies_full_clean[i]))

    c_efficiencies_peak.append(max(c_efficiencies_full[i]))
    uc_efficiencies_peak.append(max(uc_efficiencies_full[i]))
    c_efficiencies_peak_clean.append(max(c_efficiencies_full_clean[i]))
    uc_efficiencies_peak_clean.append(max(uc_efficiencies_full_clean[i]))

t_avg_peak = []
for i in range(len(c_efficiencies_peak)):
    t_avg_peak.append(i*full_day_sample_number)

#%%
# Plotting data...

# Plotting the change in heliostat reflectance over the year
plt.figure(figsize=(12, 6))
t = np.arange(num_timesteps)

for i in range(num_heliostats):
    plt.plot(t, reflectances[i, :], label=f"Heliostat {i+1}", linewidth=1.5)

plt.xlabel("Time")
plt.ylabel("Reflectance")
plt.title("Change in Heliostat Reflectances")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# Overlaying the change in reflectance with the uncorrected peak efficiencies
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t, avg_reflectances, color="blue", label="Average Reflectances")
ax1.set_xlabel("Time")
ax1.set_ylabel("Average Reflectance", color="blue")
ax1.set_title('Comparing the Average Field Reflectance to the Peak Efficiencies')
ax1.tick_params(axis="y", labelcolor="blue")

ax2 = ax1.twinx()
ax2.plot(t_avg_peak, uc_efficiencies_peak, color="red", label="Uncorrected")
ax2.plot(t_avg_peak, c_efficiencies_peak, color="green", label="Corrected")
ax2.set_ylabel("Peak Efficiency (per day)", color="black")
ax2.tick_params(axis="y", labelcolor="black")
ax2.legend(loc='lower right')
ax1.grid(True)
plt.show()

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t, avg_reflectances, color="blue", label="Average Reflectances")
ax1.set_xlabel("Time")
ax1.set_ylabel("Average Reflectance", color="blue")
ax1.set_title('Comparing the Average Field Reflectance to the Avg Efficiencies')
ax1.tick_params(axis="y", labelcolor="blue")

ax2 = ax1.twinx()  
ax2.plot(t_avg_peak, uc_efficiencies_avg, color="red", label="Uncorrected")
ax2.plot(t_avg_peak, c_efficiencies_avg, color="green", label="Corrected")
ax2.set_ylabel("Average Efficiency (per 24 hour period)", color="black")
ax2.tick_params(axis="y", labelcolor="black")
ax2.legend(loc='lower right')
ax1.grid(True)
plt.show()

# Comparing corrected clean and corrected dirty
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t_avg_peak, c_efficiencies_avg_clean, color="blue", label= "Clean Avg Efficiencies")
ax1.plot(t_avg_peak, c_efficiencies_avg, color = 'red', label = 'Soiled Avg Efficiencies')
ax1.set_xlabel("Time")
ax1.set_ylabel("Average Efficiency")
ax1.set_title('Comparing Average Field Efficiencies between Soiled and Clean Simulations')
ax1.legend()
ax1.grid(True)
plt.show()

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t_avg_peak, c_efficiencies_peak_clean, color="blue", label="Clean Peak Efficiencies")
ax1.plot(t_avg_peak, c_efficiencies_peak, color = 'red', label = 'Soiled Peak Efficiencies')
ax1.set_xlabel("Time")
ax1.set_ylabel("Peak Efficiency")
ax1.set_title('Comparing Peak Field Efficiencies between Soiled and Clean Simulations')
ax1.legend()
ax1.grid(True)
plt.show()

# Plant Efficiency Over a 24 Hour Period - Clean vs Soiled
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t[67700:67950], corrected_efficiencies_clean[67700:67950], color="red", label = 'Clean')
ax1.plot(t[67700:67950], corrected_efficiencies[67700:67950], color="green", label = 'Soiled')
ax1.set_xlabel("Time")
ax1.set_ylabel("Efficiency", color="black")
ax1.set_title('Comparing Soiled vs Clean Efficiencies throughout the Day')
ax1.tick_params(axis="y", labelcolor="black")
ax1.legend()
ax1.grid(True)
plt.show()
# %%