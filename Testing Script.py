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

wodonga_path = os.path.join(current_directory, "Wodonga Data")

# Import modules
import numpy as np
import matplotlib.pyplot as plt
import soiling_model.base_models as smb

# Define file paths
d = heliosoil_path + "/woomera_demo/"
file_params = wodonga_path + "\\parameters_wodonga_experiments.xlsx"
file_weather = wodonga_path + "\\experiment_20220220_20220226.xlsx"

# Initialize simulation data and models
sim_data = smb.simulation_inputs(file_weather, dust_type="PM10", verbose=False)

# Initialize and load site data into the physical model
physical_model = smb.physical_base()
physical_model.import_site_data_and_constants(file_params)

# Example heliostat setup. 11 Heliostats, each at a certain angle during the whole year.
num_heliostats = 11
num_timesteps = 2009
fixed_tilt_angles = np.array([5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0])
physical_model.helios.tilt = {0: np.tile(fixed_tilt_angles, (num_timesteps, 1)).T}

# Compute optical extinction weights
physical_model.helios.compute_extinction_weights(sim_data, loss_model="geometry", verbose=True)

# Calculate deposition flux to populate pdfqN
physical_model.deposition_flux(sim_data, hrz0=physical_model.hrz0, verbose=True)

# Calculate delta soiled area
sigma_dep = 0.01
physical_model.calculate_delta_soiled_area(sim_data, sigma_dep=sigma_dep, verbose=True)

# Access and analyze results
delta_soiled_area = physical_model.helios.delta_soiled_area[0]

# Generate time (hours) for the x-axis
time = np.arange(num_timesteps)

# Plot the delta_soiled_area for each heliostat
plt.figure(figsize=(12, 8))

for i in range(num_heliostats):
    plt.plot(time, delta_soiled_area[i, :], label=f"Heliostat {i+1}", linewidth=1.5)

# Add labels and title
plt.xlabel("Time (hours)")
plt.ylabel("Soiled Area (m²/m²)")
plt.title("Soiling of Heliostats Over Time")
plt.legend(loc='upper right')
plt.grid(True)

# Show the plot
plt.tight_layout()
plt.show()

