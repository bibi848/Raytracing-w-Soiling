# Get the directory of the current script
import os
import sys
current_script_path = os.path.abspath("Soiling Test Example.py") # Absolute path to a script
current_directory = os.path.dirname(current_script_path)         # Directory containing the script

# Construct the path to the HelioSoil folder and add the HelioSoil folder to the Python path
heliosoil_path = os.path.join(current_directory, "HelioSoil")
sys.path.append(heliosoil_path)
wodonga_path = os.path.join(current_directory, "Wodonga Data")
d = wodonga_path + "\\"

import pandas as pd

week1_data = pd.read_excel(d + "experiment_20220220_20220226.xlsx", sheet_name = 'Reflectance_Average')



first_header = week1_data.columns[0]
first_column_data = week1_data[first_header]

print(first_column_data)