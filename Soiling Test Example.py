
"""
Soiling Model Testing. 
This script attempts to familirise the user with the importing and manipulation of the HelioSoil soiling models, 
which later is used in conjunction with the LFR simulations to simulate the dust deposition on a CSP field.
"""

# Get the directory of the current script
import os
import sys
current_script_path = os.path.abspath(__file__)           # Absolute path to the script
current_directory = os.path.dirname(current_script_path)  # Directory containing the script

# Construct the path to the HelioSoil folder and add the HelioSoil folder to the Python path
heliosoil_path = os.path.join(current_directory, "HelioSoil")
sys.path.append(heliosoil_path)

import numpy as np
import matplotlib.pyplot as plt

# Import modules
import soiling_model.base_models as smb

# Define file paths
d = heliosoil_path + "/woomera_demo/"
file_params = d + "parameters.xlsx"
file_weather = d + "woomera_data.xlsx"

# Initialize simulation data and models
sim_data = smb.simulation_inputs(file_weather, dust_type="PM10", verbose=False)

# Initialize and load site data into the physical model
physical_model = smb.physical_base()
physical_model.import_site_data_and_constants(file_params)

# Initialize and populate the dust instance
my_dust = smb.dust()
my_dust.D[0] = np.array([1.0, 2.0, 5.0, 10.0])  # Convert to NumPy array for compatibility
my_dust.rho[0] = 2600                           # Density of dust particles in kg/m³
my_dust.pdfN[0] = np.array([100, 80, 50, 10])   # Convert to NumPy array

# Extract parameters from sim_data
wind_speed = sim_data.wind_speed[0]  # Wind speed for the first simulation file
air_temp = sim_data.air_temp[0]      # Air temperature for the first simulation file

# Call the deposition_velocity method
results = physical_model.deposition_velocity(
    dust=my_dust,
    wind_speed=wind_speed,
    air_temp=air_temp,
    hrz0=physical_model.hrz0  # Use hrz0 defined in the parameters file
)

# Generate dummy tilt data: 2 heliostats, 8760 timesteps (hourly for a year)
num_heliostats = 2
num_timesteps = 8760
physical_model.helios.tilt = {0: np.random.uniform(30, 60, (num_heliostats, num_timesteps))}

# Compute optical extinction weights
physical_model.helios.compute_extinction_weights(sim_data, loss_model="geometry", verbose=True)

# Manually set sim_data.dt if missing
if not sim_data.dt:
    sim_data.dt = {0: 3600}  # 1-hour time step

# Calculate deposition flux to populate pdfqN
physical_model.deposition_flux(sim_data, hrz0=physical_model.hrz0, verbose=True)

# Optionally define sigma_dep
sigma_dep = 0.01  # Example value

# Calculate delta soiled area
physical_model.calculate_delta_soiled_area(sim_data, sigma_dep=sigma_dep, verbose=True)

# Access and analyze results
delta_soiled_area = physical_model.helios.delta_soiled_area[0]
print(f"Delta Soiled Area (shape): {delta_soiled_area.shape}")
print(f"Delta Soiled Area (sample): {delta_soiled_area[:5, :5]}")

