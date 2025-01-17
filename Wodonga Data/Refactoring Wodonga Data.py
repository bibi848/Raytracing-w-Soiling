#%%
# Importing all Wodonga Data
# Get the directory of the current script
import os
import sys
current_script_path = os.path.abspath("Soiling Test Example.py") # Absolute path to a script
current_directory = os.path.dirname(current_script_path)         # Directory containing the script

wodonga_path = os.path.join(current_directory, "Wodonga Data")
wodonga_path += "\\"

import pandas as pd
import numpy as np

all_data = pd.read_excel("ALL_DATA_WODONGA.xlsx")

# %%
# Extracting relevant data from Wodonga dataframe
Time_data = all_data["tmsmp"][122:].reset_index(drop=True)
AirTemp_data = all_data["M1_T1"][122:].reset_index(drop=True)
WindSpeed_data = all_data["M1_WS"][122:].reset_index(drop=True)
WindDir_data = all_data["M1_WD"][122:].reset_index(drop=True)
RH_data = all_data["M1_ppRH"][122:].reset_index(drop=True)
Rain_data = all_data["M1_P"][122:].reset_index(drop=True)
DewPoint_data = all_data["S_DP"][122:].reset_index(drop=True)
PM10_data = all_data["PM10"][122:].reset_index(drop=True)
data_length = len(Time_data)

Wetbulb_data = pd.Series(np.nan, index=range(data_length))
Visibility_data = pd.Series(np.nan, index=range(data_length))
DNI_data = pd.Series(np.nan, index=range(data_length))

# %%
# Placing data into a new excel file

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

data_df = pd.DataFrame(data_dict)

data_df.to_excel("Wodonga Data Refactored.xlsx", sheet_name = 'Weather', index = False)

#%%