# Get the directory of the current script
import os
import sys
current_script_path = os.path.abspath("Soiling Test Example.py") # Absolute path to a script
current_directory = os.path.dirname(current_script_path)         # Directory containing the script

# Construct the path to the HelioSoil folder and add the HelioSoil folder to the Python path
heliosoil_path = os.path.join(current_directory, "HelioSoil")
sys.path.append(heliosoil_path)
wodonga_path = os.path.join(current_directory, "Wodonga Data")

# Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import soiling_model.base_models as smb
import soiling_model.fitting as smf
import soiling_model.utilities as smu
from matplotlib import rcParams
from copy import deepcopy
import scipy.stats as sps

# Define file paths
d = wodonga_path + "\\"
file_params = d + "parameters_wodonga_experiments.xlsx"
file_weather = d + "experiment_20220220_20220226.xlsx"

sp_save_file = f"{wodonga_path}/results/sp_fitting_results_wodonga"
cm_save_file = f"{wodonga_path}/results/cm_fitting_results_wodonga"
figure_format = ".pdf"
reflectometer_incidence_angle = 15 # angle of incidence of reflectometer
reflectometer_acceptance_angle = 12.5e-3 # half acceptance angle of reflectance measurements
second_surf = True # True if using the second-surface model. Otherwise, use first-surface
use_fitted_dust_distributions = False

train_experiments = [0] # indices for training experiments from 0 to len(files)-1
train_mirrors = ["OE_M1_T00"] # which mirrors within the experiments are used for 
k_factor = None # None sets equal to 1.0, "import" imports from the file
dust_type = "PM10"
epsilon = 1e-3
M = 100 # number of Monte Carlo simulations for deposition rate histogram

files,all_intervals,exp_mirrors,all_mirrors = smu.get_training_data(d,"experiment_")
orientation = [ [s[1] for s in mirrors] for mirrors in exp_mirrors]

# Feb 2022 (first experiment --- remove last three days after rain started)
all_intervals[0][0] = np.datetime64('2022-02-20T16:20:00')
all_intervals[0][1] = np.datetime64('2022-02-23T17:40:00')

# April 2022 (remove nothing, first/last measurement)
all_intervals[1][0] = np.datetime64('2022-04-21T11:00:00')
all_intervals[1][1] = np.datetime64('2022-04-27T08:30:00')

# Feb 2023 (most recent experiment --- remove very dirty days)
all_intervals[2][0] = np.datetime64('2023-02-09T15:00:00')
all_intervals[2][1] = np.datetime64('2023-02-14T09:45:00')

testing_intervals = all_intervals
        
Nfiles = len(files)
extract = lambda x,ind: [x[ii] for ii in ind]
files_train = extract(files,train_experiments)
training_intervals = extract(all_intervals,train_experiments)
testing_intervals = list(all_intervals)
t = [t for t in train_experiments]
plot_title = "Training: "+str(train_mirrors)+", Exp: "+str(t)

# Import and Plot Training Data
imodel = smf.semi_physical(file_params)
imodel_constant = smf.constant_mean_deposition(file_params)
sim_data_train = smb.simulation_inputs( files_train,
                                        k_factors=k_factor,
                                        dust_type=dust_type
                                        )
reflect_data_train = smb.reflectance_measurements(  files_train,
                                                    sim_data_train.time,
                                                    number_of_measurements=9.0,
                                                    reflectometer_incidence_angle=reflectometer_incidence_angle,
                                                    reflectometer_acceptance_angle=reflectometer_acceptance_angle,
                                                    import_tilts=True,
                                                    column_names_to_import=train_mirrors
                                                    )

# Trim data and plot
sim_data_train,reflect_data_train = smu.trim_experiment_data(   sim_data_train,
                                                                reflect_data_train,
                                                                training_intervals 
                                                            )
                                                            
sim_data_train,reflect_data_train = smu.trim_experiment_data(   sim_data_train,
                                                                reflect_data_train,
                                                                "reflectance_data" 
                                                            )
for ii,experiment in enumerate(train_experiments):
    fig,ax = smu.plot_experiment_data(sim_data_train,reflect_data_train,ii)
    fig.suptitle(f"Training Data for file {files[experiment]}")

# Set mirror angles and get extinction weights
imodel.helios_angles(sim_data_train,reflect_data_train,second_surface=second_surf)
imodel.helios.compute_extinction_weights(sim_data_train,imodel.loss_model,
                                         verbose=False,options={'grid_size_x':1000})
fig_weights,ax_weights = imodel.helios.plot_extinction_weights(sim_data_train,fig_kwargs={'figsize':(5,7)})
ext_weights = imodel.helios.extinction_weighting[0].copy()

imodel_constant.helios_angles(sim_data_train,reflect_data_train,second_surface=second_surf)
file_inds = np.arange(len(files_train))
imodel_constant = smu.set_extinction_coefficients(imodel_constant,ext_weights,file_inds)

# Fit semi-physical model & plot on training data
log_param_hat,log_param_cov = imodel.fit_mle(   sim_data_train,
                                        reflect_data_train,
                                        transform_to_original_scale=False)

s = np.sqrt(np.diag(log_param_cov))
param_ci = log_param_hat + 1.96*s*np.array([[-1],[1]])
lower_ci = imodel.transform_scale(param_ci[0,:])
upper_ci = imodel.transform_scale(param_ci[1,:])
param_hat = imodel.transform_scale(log_param_hat)
hrz0_mle,sigma_dep_mle = param_hat
print(f'hrz0: {hrz0_mle:.2e} [{lower_ci[0]:.2e}, {upper_ci[0]:.2e}]')
print(f'\sigma_dep: {sigma_dep_mle:.2e} [{lower_ci[1]:.2e},{upper_ci[1]:.2e}] [p.p./day]')

imodel.update_model_parameters(param_hat)
imodel.save(sp_save_file,
            log_p_hat=log_param_hat,
            log_p_hat_cov=log_param_cov,
            training_simulation_data=sim_data_train,
            training_reflectance_data=reflect_data_train)

_,_,_ = imodel.plot_soiling_factor( sim_data_train,
                            reflectance_data=reflect_data_train,
                            figsize=(10,10),
                            reflectance_std='measurements',
                            save_path=f"{wodonga_path}/results/wodonga_semi_physical_training",
                            fig_title="On Training Data (semi-physical)",
                            orientation_strings=orientation    )

# Load and trim total data
sim_data_total = smb.simulation_inputs( files,
                                        k_factors=k_factor,
                                        dust_type=dust_type
                                        )

reflect_data_total = smb.reflectance_measurements(  files,
                                                    sim_data_total.time,
                                                    number_of_measurements=9.0,
                                                    reflectometer_incidence_angle=reflectometer_incidence_angle,
                                                    reflectometer_acceptance_angle=reflectometer_acceptance_angle,
                                                    import_tilts=True,
                                                    column_names_to_import=None
                                                    )
sim_data_total,reflect_data_total = smu.trim_experiment_data(   sim_data_total,
                                                                reflect_data_total,
                                                                testing_intervals 
                                                            )

sim_data_total,reflect_data_total = smu.trim_experiment_data(   sim_data_total,
                                                                reflect_data_total,
                                                                "reflectance_data" 
                                                            )

sim_data_total,reflect_data_total = smu.trim_experiment_data(   sim_data_total,
                                                                reflect_data_total,
                                                                "simulation_inputs" 
                                                            )

# Plot Experiments
for ii,experiment in enumerate(sim_data_total.dt.keys()):
    fig,ax = smu.plot_experiment_data(sim_data_total,reflect_data_total,ii)
    fig.suptitle(f"Testing Data for file {files[experiment]}")

