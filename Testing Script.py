
# Get the directory of the current script
import os
import sys
current_script_path = os.path.abspath("Soiling Test Example.py") # Absolute path to a script
current_directory = os.path.dirname(current_script_path)         # Directory containing the script

# Construct the path to the HelioSoil folder and add the HelioSoil folder to the Python path
heliosoil_path = os.path.join(current_directory, "HelioSoil")
sys.path.append(heliosoil_path)
wodonga_path = os.path.join(current_directory, "Wodonga Data")
wodonga_path += "\\"

import pandas as pd

all_data = pd.read_excel(wodonga_path + "ALL_DATA_WODONGA.xlsx")


