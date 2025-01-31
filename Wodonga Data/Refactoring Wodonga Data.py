#%%
# Importing all Wodonga Data
"""
This script takes ALL_DATA_WODONGA.xlsx and extracts the data required to be placed in a new excel file which can be 
read by the HelioSoil scripts. 
"""
# Get the directory of the current script
import os
current_script_path = os.path.abspath(__file__)          # Absolute path to a script
wodonga_directory = os.path.dirname(current_script_path) # Directory containing the script
HelioSoil_path = wodonga_directory.replace(r"\Wodonga Data", "\\HelioSoil\\woomera_demo\\woomera_data.xlsx")

import pandas as pd
import numpy as np
from datetime import datetime
import datetime as dt
import time
import copy

all_data_t1 = pd.read_excel(wodonga_directory + "\\ALL_DATA_WODONGA.xlsx")
# all_data = pd.read_excel(wodonga_directory + "\\test set.xlsx")

# Only taking data we need to make computation faster
Time_data = all_data_t1["tmsmp"].reset_index(drop=True)
AirTemp_data = all_data_t1["M1_T1"].reset_index(drop=True)
WindSpeed_data = all_data_t1["M1_WS"].reset_index(drop=True)
WindDir_data = all_data_t1["M1_WD"].reset_index(drop=True)
RH_data = all_data_t1["M1_ppRH"].reset_index(drop=True)
Rain_data = all_data_t1["M1_P"].reset_index(drop=True)
DewPoint_data = all_data_t1["S_DP"].reset_index(drop=True)
PM10_data = all_data_t1["PM10"].reset_index(drop=True)
data_length = len(Time_data)

Wetbulb_data = pd.Series(0, index=range(data_length))
Visibility_data = pd.Series(0, index=range(data_length))
DNI_data = pd.Series(0, index=range(data_length))

data_dict = {
    "Time" : Time_data,
    "AirTemp" : AirTemp_data,
    "WindSpeed" : WindSpeed_data,
    "WindDir" : WindDir_data,
    "RH" : RH_data,
    "Rain" : Rain_data,
    "DewPoint" : DewPoint_data,
    "WetBulb" : Wetbulb_data,
    "Visibility" : Visibility_data,
    "PM10" : PM10_data,
    "DNI" : DNI_data
}

all_data_t2 = pd.DataFrame(data_dict)
#%%
# Functions for data refactoring...

def linear_interpolation(x0, y0, x1, y1, x):
    return (y0*(x1 - x) + y1*(x - x0)) / (x1 - x0)

def fill_gap(index, all_data, time_diff, middle_time):
    interpolated_results = [middle_time]
    for column in all_data.columns[1:]:
        interpolated_results.append(linear_interpolation(0, all_data[column][index], time_diff, all_data[column][index+1], time_diff/2))
    return interpolated_results

def round_timestamp(timestamp):
    rounded_timestamp = timestamp - pd.Timedelta(minutes=timestamp.minute % 5, seconds=timestamp.second, microseconds=timestamp.microsecond)
    if timestamp.minute % 5 >= 3:
        rounded_timestamp += pd.Timedelta(minutes=5)
    return rounded_timestamp

def check_5_minute_interval(timestamp, counter):
    if timestamp.minute % 5 != 0 or timestamp.second != 0:
        counter += 1
        return round_timestamp(timestamp), counter
    else:
        return timestamp, counter

#%%
# Filling 2 month gap with previous data
# No data is available between 2022-05-20 14:45:00 (103081) and 2022-08-09 10:05:46 (103082).
# Will be replaced by 2021-05-20 14:45:00 (11) and 2021-08-09 10:05:00 (21583)

filler_data = copy.deepcopy(all_data_t2.iloc[12:21583])
filler_data.loc[:, 'Time'] = filler_data['Time'].apply(lambda x: x.replace(year=2022))

all_data_part1 = all_data_t2.iloc[:103082]
all_data_part2 = all_data_t2.iloc[103082:]
all_data = pd.concat([all_data_part1, filler_data, all_data_part2], ignore_index = True)
print('Two month gap has been filled')

#%%
# Remove all data to first instance of midnight
for i in range(len(all_data["Time"])):
    tmsmp_value = all_data["Time"][i]
    if isinstance(tmsmp_value, (pd.Timestamp, dt.datetime)):
        if tmsmp_value.time() == dt.time(0, 0):
            break
all_data = all_data.iloc[i:].reset_index(drop = True)
print(i, 'rows have been removed off the start')

#%%
# Ensuring all data is in 5 minute intervals
incorrect_times = 0
for i in range(len(all_data["Time"])):
    all_data.loc[i, 'Time'], incorrect_times = check_5_minute_interval(all_data["Time"][i], incorrect_times)

print(incorrect_times, 'times were not in 5 minute intervals')

#%%
# Filling the gaps
gaps_filled = 0
i = 0
while i < len(all_data["Time"]):
    if i != len(all_data["Time"])-1:

        current_time = all_data["Time"][i]
        next_time = all_data["Time"][i+1]
        time_diff = (next_time - current_time)

        if time_diff.total_seconds() != 300.0:
            print(current_time)

            middle_time = round_timestamp(current_time + (next_time - current_time) / 2)

            gaps_filled += 1
            interpolated_results = fill_gap(i, all_data, time_diff.total_seconds(), middle_time)

            new_data = pd.DataFrame([interpolated_results], columns = all_data.columns)

            all_data_part1 = all_data.iloc[:i + 1]
            all_data_part2 = all_data.iloc[i + 1:]

            all_data = pd.concat([all_data_part1, new_data, all_data_part2], ignore_index = True)
            i -= 1
            
    i += 1

print(gaps_filled,'gaps were filled by interpolation')

# %%
# Getting Dust sheet from Woomera Data
woomera_data = pd.read_excel(HelioSoil_path, sheet_name= "Dust")

# %%
# Placing data into a new excel file

with pd.ExcelWriter(wodonga_directory + "\\Wodonga Data Refactored.xlsx", engine='openpyxl') as writer:
    woomera_data.to_excel(writer, sheet_name='Dust', index=False)
    all_data.to_excel(writer, sheet_name='Weather', index=False)

#%%