# Wodonga Data

A folder within the Raytracing-w-Soiling repository which analyses the effect of soiling on the effectiveness of a Linear Fresnel Reflector (LFR) plant in Wodonga, VIC Australia.

## Summary
This folder provides simulations and tools to analyse the effect of the soiling models developed in the HelioSoil [1] repository on the effectiveness of a LFR plant in Wodonga simulated using PySolTrace by NREL [2]. The data for the case study has been collected by the Queensland University of Technology.

All scripts contain comments to explain their use, but a brief explanation for each is below. 

Folders:
* Training Data: This contains the training data for the hrz0 fitting model from HelioSoil.
<br><br>
Files:
* `ALL_DATA_WODONGA.xlsx`: Collected data from the Wodonga site provided by the Queensland University of Technology.
`README.md`: The README file you are currently reading.
* `Refactoring Wodonga Data.py`: The `ALL_DATA_WODONGA.xlsx` does not present the data in the correct manner. This Python script therefore extracts the required data and also trims the initial data so that the data starts at midnight. 
* `Wodonga Data Refactored.xlsx`: This is the script produced from `Refactoring Wodonga Data.py`. It contains the weather and dust data which is then used in `Wodonga Plant Soiling Analysis.py`.
* `Wodonga Graphing.py`: This script is used to analyse and graph all data produced from the simulations.
* `Wodonga Plant Soiling Analysis.py`: This script computes the soiling analysis for the designed LFR plant. It produces `Wodonga Soiled Data.xlsx`, but this file is too large to currently be uploaded to GitHub. To access it, use this repository and it will be created on your local system.
* `Wodonga Ray Tracing.py`: From the the soiling analysis, the LFR plant is subjected to the ray tracing simulations. This produces `Wodonga Raytrace Results.csv`.
* `Wodonga Simulation Parameters.csv`: This contains the physical parameters for the design of the LFR plant.
* `parameters_wodonga_experiements.xlsx`: This contains other parameters relating to the location of the Wodonga plant, such as some heliostat and air constants.

<br><br>
 ## Main Workflow Structure
 After setting the physical plant parameters in the `Wodonga Simulation Parameters.csv` file, the `Refactoring Wodonga Data.py` script is run to ensure that the data used from `ALL_DATA_WODONGA.xlsx` is in the correct format. From there, the `Wodonga Plant Soiling Analysis.py` script is run to compute all soiling parameters from the data set on the LFR plant. This data is passed to `Wodonga Ray Tracing.py` to get a measure of the plant's performance which is visualised and analysed in `Wodonga Graphing.py`. 
 
 


