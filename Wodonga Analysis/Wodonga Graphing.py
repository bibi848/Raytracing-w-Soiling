#%% 
# Importing Modules and Data...

"""
This script takes the data produced both by the ray tracing and soiling analysis scripts and provides a space where the data can
be analysed and visulised. 
"""

# Finding the script's directory to then find the data csv files produced in previous scripts.
import os
import sys
current_script_path = os.path.abspath(__file__)
wodonga_path = os.path.dirname(current_script_path)
current_directory = wodonga_path.replace(r"\Wodonga Analysis", "")
sys.path.append(current_directory)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, time
from raytracing_soiling_functions import import_simulation_parameters

csv_path = os.path.join(current_directory, "Wodonga Analysis\\Wodonga Data\\Wodonga Soiled Data.csv")
df_soiled_data = pd.read_csv(csv_path)
csv_path = os.path.join(current_directory, "Wodonga Analysis\\Wodonga Data\\Wodonga Clean Data.csv")
df_clean_data = pd.read_csv(csv_path)
csv_path = os.path.join(current_directory, "Wodonga Analysis\\Wodonga Data\\Wodonga Raytrace Results.csv")
df_raytrace_results = pd.read_csv(csv_path)
csv_path = os.path.join(current_directory, "Wodonga Analysis\\Wodonga Data\\Wodonga Raytrace Results - Clean.csv")
df_raytrace_results_clean = pd.read_csv(csv_path)

# Extracting the data from the csv file created in Wodonga Plant Soiling Analysis.py
azimuths_deg   = df_soiled_data['Azimuth [deg]'].to_numpy()
azimuths_rad   = np.deg2rad(azimuths_deg)
elevations_deg = df_soiled_data['Elevation [deg]'].to_numpy()
elevations_rad = np.deg2rad(elevations_deg)
num_timesteps  = len(df_soiled_data['Date'].to_numpy())
DNI_day_sums   = df_soiled_data["DNI Sum Per Day"].to_numpy()
DNI_all        = df_soiled_data["DNI [W/m2]"].to_numpy()
energy_demand  = df_soiled_data["Energy Demand [kWh]"].to_numpy()

csv_path = current_directory + "\\Wodonga Analysis\\Wodonga Data\\Wodonga Simulation Parameters.csv"
lat, lon, hour_offset, receiver_height, receiver_length, receiver_diameter, panel_length, panel_width, panel_height, panel_spacing, number_of_modules, panels_per_module, slope_error, specularity_error = import_simulation_parameters(pd.read_csv(csv_path))
num_heliostats = number_of_modules * panels_per_module
aperture = num_heliostats * panel_length * panel_width

tilt_header_list = []
cleanliness_header_list = []
cumulative_soiling_header_list = []
for i in range(num_heliostats):
    tilt_header = f"Heliostat tilt [deg] {i+1}"
    cleanliness_header = f"Heliostat Cleanliness {i+1}"
    cumulative_soiling_header = f"Cumulative Soiled Area [m2/m2] {i+1}"

    tilt_header_list.append(tilt_header)
    cleanliness_header_list.append(cleanliness_header)
    cumulative_soiling_header_list.append(cumulative_soiling_header)

tilts_deg = df_soiled_data[tilt_header_list].to_numpy().T
cleanlinesses = df_soiled_data[cleanliness_header_list].to_numpy().T
cumulative_soiled_area = df_soiled_data[cumulative_soiling_header_list].to_numpy().T

# Extracting the data from the csv file created in Wodonga Ray Tracing.py
optical_efficiencies = df_raytrace_results["Optical efficiency"].to_numpy()
field_efficiencies = df_raytrace_results["Field efficiency"].to_numpy()
DNI_efficiencies = df_raytrace_results["DNI corrected efficiency"].to_numpy()
optical_efficiencies_clean = df_raytrace_results_clean["Optical efficiency"].to_numpy()
field_efficiencies_clean = df_raytrace_results_clean["Field efficiency"].to_numpy()
DNI_efficiencies_clean = df_raytrace_results_clean["DNI corrected efficiency"].to_numpy()

#%% 
# Manipulating data...

# Average Field Reflectance
avg_cleanliness = []
for i in range(num_timesteps):
    avg = []
    for p in range(num_heliostats):
        avg.append(cleanlinesses[p][i])
    avg_cleanliness.append(sum(avg) / len(avg))

o_efficiencies_full = []
f_efficiencies_full = []
D_efficiencies_full = []
o_efficiencies_full_clean = []
f_efficiencies_full_clean = []
D_efficiencies_full_clean = []

o_efficiencies_day = []
f_efficiencies_day = []
D_efficiencies_day = []
o_efficiencies_day_clean = []
f_efficiencies_day_clean = []
D_efficiencies_day_clean = []

for i in range(len(optical_efficiencies)):

    o_efficiency = optical_efficiencies[i]
    f_efficiency = field_efficiencies[i]
    D_efficiency = DNI_efficiencies[i]
    o_efficiency_clean = optical_efficiencies_clean[i]
    f_efficiency_clean = field_efficiencies_clean[i]
    D_efficiency_clean = DNI_efficiencies_clean[i]

    if o_efficiency != 0:

        o_efficiencies_day.append(o_efficiency)
        f_efficiencies_day.append(f_efficiency)
        D_efficiencies_day.append(D_efficiency)
        o_efficiencies_day_clean.append(o_efficiency_clean)
        f_efficiencies_day_clean.append(f_efficiency_clean)
        D_efficiencies_day_clean.append(D_efficiency_clean)

        if i < len(optical_efficiencies)-1 and optical_efficiencies[i+1] < 0.01:

            o_efficiencies_full.append(o_efficiencies_day)
            f_efficiencies_full.append(f_efficiencies_day)
            D_efficiencies_full.append(D_efficiencies_day)
            o_efficiencies_full_clean.append(o_efficiencies_day_clean)
            f_efficiencies_full_clean.append(f_efficiencies_day_clean)
            D_efficiencies_full_clean.append(D_efficiencies_day_clean)

            o_efficiencies_day = []
            f_efficiencies_day = []
            D_efficiencies_day = []
            o_efficiencies_day_clean = []
            f_efficiencies_day_clean = []
            D_efficiencies_day_clean = []

o_efficiencies_avg = []
f_efficiencies_avg = []
D_efficiencies_avg = []
o_efficiencies_avg_clean = []
f_efficiencies_avg_clean = []
D_efficiencies_avg_clean = []

o_efficiencies_peak = []
f_efficiencies_peak = []
D_efficiencies_peak = []
o_efficiencies_peak_clean = []
f_efficiencies_peak_clean = []
D_efficiencies_peak_clean = []

for i in range(len(o_efficiencies_full)):
    o_efficiencies_avg.append(np.mean(o_efficiencies_full[i]))
    f_efficiencies_avg.append(np.mean(f_efficiencies_full[i]))
    D_efficiencies_avg.append(sum(D_efficiencies_full[i]) / DNI_day_sums[i])
    o_efficiencies_avg_clean.append(np.mean(o_efficiencies_full_clean[i]))
    f_efficiencies_avg_clean.append(np.mean(f_efficiencies_full_clean[i]))
    D_efficiencies_avg_clean.append(sum(D_efficiencies_full_clean[i]) / DNI_day_sums[i])

    o_efficiencies_peak.append(max(o_efficiencies_full[i]))
    f_efficiencies_peak.append(max(f_efficiencies_full[i]))
    D_efficiencies_peak.append(max(D_efficiencies_full[i]))
    o_efficiencies_peak_clean.append(max(o_efficiencies_full_clean[i]))
    f_efficiencies_peak_clean.append(max(f_efficiencies_full_clean[i]))

# As the efficiency has been calculated for every 5 minutes, it is broken down into blocks of 288, which represent 24 hour periods.
full_day_sample_number = (24 * 60) / 5  # Number of minutes in a full 24 hour period, divided by the sample time step (5 minutes)

t_avg_peak = []
for i in range(len(o_efficiencies_peak)):
    t_avg_peak.append(i*full_day_sample_number)
for i in range(len(t_avg_peak)):
    t_avg_peak[i] = t_avg_peak[i] * (5/60)

#%%
# Energy Analysis: Energy Generated

# Finding the indexes of the first of each month
month_indexes = []
scale_factor = 30  # If simulated field is smaller than actual field required, this is the multiplication factor.

for i in range(num_timesteps):
    date_obj = datetime.fromisoformat(df_soiled_data["Date"][i])
    day = date_obj.day
    is_midnight = date_obj.time()

    if day == 1 and is_midnight == time(0,0,0):
        month_indexes.append(i)

# Finding the energy generated each timestep
def calc_power(field_efficiency, DNI, field_area, receiver_efficiency):
    return field_efficiency * DNI * field_area * receiver_efficiency

receiver_efficiency = 0.9
energy_generated = []
energy_generated_clean = []
power_generated = []

for i in range(num_timesteps):

    DNI_ts = DNI_all[i]
    field_efficiency_ts = field_efficiencies[i]
    field_efficiency_ts_clean = field_efficiencies_clean[i]

    power = (scale_factor * calc_power(field_efficiency_ts, DNI_ts, aperture, receiver_efficiency)) / 1000 # kW
    power_clean = (scale_factor * calc_power(field_efficiency_ts_clean, DNI_ts, aperture, receiver_efficiency)) / 1000

    ener = power / 12 # kWh over a period of 5 minutes
    ener_clean = power_clean / 12

    energy_generated.append(ener)
    energy_generated_clean.append(ener_clean)
    power_generated.append(power)

# Grouping the energy generated per month
energy_generated_months = []
energy_generated_clean_months = []

for i in range(len(month_indexes)-1):
    start = month_indexes[i]
    end = month_indexes[i+1]

    month_sum = sum(energy_generated[start:end])
    month_sum_clean = sum(energy_generated_clean[start:end])

    energy_generated_months.append(month_sum)
    energy_generated_clean_months.append(month_sum_clean)


labels = ["Jun '21", "Jul '21", "Aug '21", "Sep '21", "Oct '21", "Nov '21", "Dec '21",
          "Jan '22", "Feb '22", "Mar '22", "Apr '22", "May '22", "Jun '22", "Jul '22",
          "Aug '22", "Sep '22", "Oct '22", "Nov '22", "Dec '22", "Jan '23", 
          ]

#%%
# Plotting data...

t = np.arange(num_timesteps)
for i in range(len(t)):
    t[i] = t[i] * (5/60)

# The change in cumulative soiled area for each heliostat over the time period.

fig, ax1 = plt.subplots(figsize=(12,6))

for i in range(num_heliostats):
    ax1.plot(t, cumulative_soiled_area[i, :], linewidth=1.5)

ax1.set_xlabel("Time [Hours]")
ax1.set_ylabel("Cumulative Soiled Area [m2/m2]")
ax1.set_title("Change in Panel Cumulative Soiled Area")
ax1.grid(True)
plt.show()

# The change in cleanliness for each heliostat over the time period.
fig, ax1 = plt.subplots(figsize=(12, 6))

for i in range(num_heliostats):
    ax1.plot(t, cleanlinesses[i, :], linewidth=1.5)

ax1.set_xlabel("Time [Hours]")
ax1.set_ylabel("Cleanliness")
ax1.set_title("Change in Panel Cleanliness")
ax1.grid(True)
plt.show()

# Overlaying the average field cleanliness with the average optical and field efficiencies.
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t, avg_cleanliness, color="grey")
ax1.set_xlabel("Time [Hours]", fontsize = 16)
ax1.set_ylabel("Average Cleanliness", color="black", fontsize = 16)
#ax1.set_title('Comparing the Average Field Cleanliness to the Average Optical and Field Efficiencies')
ax1.tick_params(axis="y", labelcolor="black")
ax2 = ax1.twinx()
#ax2.plot(t_avg_peak, o_efficiencies_avg, color="red", label="Optical")
ax2.plot(t_avg_peak, f_efficiencies_avg, color="red")
ax2.set_ylabel("Daily Average Efficiency", color="red", fontsize = 16)
ax2.tick_params(axis="y", labelcolor="red")
ax1.grid(True)
plt.show()

# Comparing the clean and dirty average field efficiencies
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t_avg_peak, f_efficiencies_avg_clean, color= 'blue', label= 'Clean')
ax1.plot(t_avg_peak, f_efficiencies_avg, color = 'red', label = 'Soiled')
ax1.axvline(8500, color = 'red', linestyle = ':', linewidth = 2)
ax1.axvline(10500, color = 'red', linestyle = ':', linewidth = 2)
ax1.set_xlabel("Time [Hours]")
ax1.set_ylabel("Daily Average Efficiency")
ax1.set_title('Comparing the Clean and Dirty Average Field Efficiencies')
ax1.legend()
ax1.grid(True)
plt.show()

# Comparing the clean and dirty average DNI corrected efficiencies
fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t_avg_peak, D_efficiencies_avg_clean, color = 'blue', label = 'Clean')
ax1.plot(t_avg_peak, D_efficiencies_avg, color = 'red', label = 'Soiled')
ax1.axvline(8500, color = 'red', linestyle = ':', linewidth = 2)
ax1.axvline(10500, color = 'red', linestyle = ':', linewidth = 2)
ax1.set_xlabel("Time [Hours]")
ax1.set_ylabel("Daily Average Efficiency")
ax1.set_title('Comparing the Clean and Dirty Average DNI corrected Efficiencies')
ax1.legend()
ax1.grid(True)
plt.show()

# %%
# DNI across the timeframe.
csv_path = os.path.join(current_directory, "Wodonga Analysis\\Wodonga Data\\Wodonga Soiled Data.csv")
df = pd.read_csv(csv_path)

DNI = df["DNI [W/m2]"].to_numpy()

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t, DNI, color="blue")
ax1.set_xlabel("Time [Hours]")
ax1.set_ylabel("DNI [W/m2]", color="blue")
ax1.set_title('DNI Variation')
ax1.tick_params(axis="y", labelcolor="blue")
ax1.grid(True)
plt.show()


# %%
# Energy Generation
fig, ax1 = plt.subplots(figsize=(12, 5))
x = np.arange(len(labels))

ax1.bar(x, energy_generated_clean_months, color='gray', width=0.4, alpha=0.5, label="Potential")
ax1.bar(x, energy_generated_months, color='blue', width=0.4, label="Actual")

ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=45)
ax1.set_ylabel("Total Energy per Month [kWh]")
ax1.set_title("Energy Production Across Different Months")
ax1.legend()

plt.show()

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t, energy_generated_clean, color= 'blue', label= 'Clean')
ax1.plot(t, energy_generated, color = 'red', label = 'Soiled')
ax1.set_xlabel("Time [Hours]")
ax1.set_ylabel("Energy Generated [kWh]")
ax1.set_title('Comparing the Clean and Dirty Energy Generations')
ax1.legend()
ax1.grid(True)
plt.show()

#%%
# Energy Usage

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t, energy_demand, color= 'grey', label = 'Demand')
ax1.plot(t, energy_generated_clean, color = 'blue', label = '"Clean" Supply')
ax1.plot(t, energy_generated, color= 'red', label = '"Soiled" Supply')
ax1.set_xlabel("Time [Hours]", fontsize = 16)
ax1.set_ylabel("Energy [kWh]", fontsize = 16)
# ax1.set_title('Energy Production vs Energy Demand')
ax1.legend(fontsize = 12)
ax1.set_ylim(500, 2300)
ax1.grid(True)
plt.show()

# %%
# Thermal Storage
storage = []
count = 0

for i in range(len(energy_generated)):
    count += energy_generated[i]
    count -= energy_demand[i]
    storage.append(count)

fig, ax1 = plt.subplots(figsize=(12, 6))
ax1.plot(t[:150], energy_generated[:150], color= 'red', label = 'Generated')
ax1.plot(t[:150], energy_demand[:150], color= 'blue', label = 'Demand')

ax2 = ax1.twinx()
ax2.plot(t[:150], storage[:150], color= 'grey', label = 'Storage')

ax1.set_xlabel("Time [Hours]", fontsize = 16)
ax1.set_ylabel("", fontsize = 16)
# ax1.set_title('Energy Production vs Energy Demand')
ax1.legend(fontsize = 12)
# ax1.set_ylim(500, 2300)
ax1.grid(True)
plt.show()