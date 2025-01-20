#%%
# Importing all Wodonga Data
# Get the directory of the current script
import os
import sys
current_script_path = os.path.abspath(__file__)          # Absolute path to a script
wodonga_directory = os.path.dirname(current_script_path) # Directory containing the script
HelioSoil_path = wodonga_directory.replace(r"\Wodonga Data", "\\HelioSoil\\woomera_demo\\woomera_data.xlsx")

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
# Getting Dust sheet from Woomera Data
woomera_data = pd.read_excel(HelioSoil_path, sheet_name= "Dust")

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

with pd.ExcelWriter("Wodonga Data Refactored.xlsx", engine='openpyxl') as writer:
    woomera_data.to_excel(writer, sheet_name='Dust', index=False)
    data_df.to_excel(writer, sheet_name='Weather', index=False)

#%%