
"""
Soiling Model Testing. 
This script attempts to familirise the user with the importing and manipulation of the HelioSoil soiling models, 
which later is used in conjunction with the LFR simulations to simulate the dust deposition on a CSP field.
"""

import os
import sys

# Get the directory of the current script
current_script_path = os.path.abspath(__file__)           # Absolute path to the script
current_directory = os.path.dirname(current_script_path)  # Directory containing the script

# Construct the path to the HelioSoil folder and add the HelioSoil folder to the Python path
heliosoil_path = os.path.join(current_directory, "HelioSoil")
sys.path.append(heliosoil_path)
print(heliosoil_path)

# Import models
import soiling_model.base_models as smb
import soiling_model.field_models as smf
import soiling_model.utilities as smu
from soiling_model.base_models import simulation_inputs

# Parameters used for the model
d = heliosoil_path + "/woomera_demo/"
file_params = d+"parameters.xlsx"
file_weather = d+'woomera_data.xlsx'

file_SF = d+'SF_woomera_SolarPILOT.csv'             # solar field of 48 sectors located in Woomera
climate_file = d+'woomera_location_modified.epw'    # only used for optical efficiency computation

n_az, n_rad = (6,8)

sim_data = smb.simulation_inputs(file_weather,dust_type="PM10", verbose = False)
imodel = smf.field_model(file_params, file_SF, num_sectors = (n_az, n_rad))
plant = smf.central_tower_plant()
plant.import_plant(file_params)
print(imodel)
# imodel.helios.sector_plot()

# imodel.sun_angles(sim_data)
# imodel.helios_angles(plant)

# imodel.compute_acceptance_angles(plant)    
# imodel.helios.compute_extinction_weights(sim_data,imodel.loss_model,verbose=True,options={'grid_size_x':250})

# imodel.deposition_flux(sim_data)
# imodel.adhesion_removal(sim_data)
# imodel.calculate_delta_soiled_area(sim_data)

# airT = 20
# windS = 2.0
# experiment = 0
# heliostat_id = 0
# imodel.plot_area_flux(sim_data,experiment,heliostat_id,airT,windS,tilt=0.0)
