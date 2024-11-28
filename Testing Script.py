
# Commercial Trapezoidal 4 Stages Lookup Table Soiled, 14 Day Cleaning Rotation

# System Setup
# Coordinate System, XYZ: X = West, Y = Perpendicular to the ground, Z = North.
# The 1st mirror is the East-most mirror, and the 11th mirror is the West-most mirror.

# FUNCTIONS

import numpy as np
import pysolar.solar as solar
import pysolar.radiation as radiation
from pysoltrace import PySolTrace, Point
from optical_geometrical_setup import stg_surface
import datetime as dt
import pandas as pd
from datetime import datetime, timedelta
import time
import os

timer_start = time.time()
def panel_normal(tilt_rad): # returns the normal vector of a panel with a certain tilt
    return np.array([np.sin(tilt_rad), np.cos(tilt_rad), 0])

# File Location for Saving
folder = os.path.dirname(os.path.abspath(__file__)) + "/CSV Result Files"
file_name = 'Commerical Tr 4 Stages Example.csv'
file_path = os.path.join(folder, file_name)

# Location: Woomera
lat, lon = -31.2,136.816667
timezoneOffset = dt.timedelta(hours = 9.5)

# Time Duration: Measuring across a whole year
# It starts on the 1st Jan 2023, and ends on the 1st Jan 2024, with a date for every hour of the year.
start_date = datetime(2023, 1, 1,hour=6,minute=0,second=0,tzinfo=dt.timezone(timezoneOffset)) # Define the start date
end_date = start_date + timedelta(days=9) # Calculate the end date 
datetime_list = []                          # Create a list of datetime objects
current_date = start_date                   # Loop to generate datetime objects from start_date to end_date
while current_date <= end_date:
    datetime_list.append(current_date)
    current_date += timedelta(hours=1)

# Simulation Parameters
number_hits = 1e3  # 1 million rays 
sunshape_flag = False
sfcerr_flag = False
cleaning_frequency = 3 # number of days between each clean

# Geometry of Panels
total_width = 7.5    # [m]
spacing = 0.2        # [m]
panel_width = 0.5    # [m]
panel_length = 12.18 # [m]
mirrors_height = 0.5 # [m]
panel_positions = np.arange(-3.5, 3.75, panel_width + spacing) # Question?
N_panels = len(panel_positions)                                 # Number of panels, 11
reflectivity = [0.95] * N_panels
reflectivity_rec = 0.06 # Absorption of 94%.
slope_error = 0.1       # [mrad]
specularity_error = 0.1 # [mrad]

# Receiver Geometry
receiver_height = 4.5          # Top of cylinder from the ground 
receiver_diameter = 0.15       # [m]
receiver_length = panel_length # [m]
receiver_optics = "Tube"       # Must match the SolTrace Optics Tab.

# All results are kept in a DataFrame
# Headers are shown below, and then panel characteristics are added (1 to 11)
lookup_table_field = ["Date", "Azimuth", "Elevation", "Optical Efficiency", "IAM", "Day/Night", 
                      "Month", "Day", "Hour", "rho_avg", "DNI", "E_sun", "E_sun_real", "E_rec"]

panel_characteristics = ['theta_', 'tilt_angle_', 'rho_', 'a_loss_', 'g_factor_']
for characteristic in panel_characteristics:
    for num in range(len(panel_positions)):
        lookup_table_field.append(characteristic + str(num+1))

# The lookup table is then initialised into a data frame, with all the panels initially having a 0 
# absorbtion loss. 
lookup_table_soiled = pd.DataFrame(columns = lookup_table_field)
lookup_table_loc = 0
a_loss = np.zeros(len(panel_positions))

# Incidence angle Evaluation
towerToMirror = np.array([np.zeros(len(panel_positions)),panel_positions,                    # Relative postion of tower from mirror.
                          np.ones((len(panel_positions)))*(receiver_height-mirrors_height)]) # Left-handed reference system.
towerToMirror_norm = towerToMirror / np.sqrt(np.sum(towerToMirror**2,axis=0))                # Normalised distance of the tower to mirror.

# Stage 1
stg1_height = receiver_height + 0.1
stg1_length = panel_length
stg1_width = panel_width * N_panels + spacing * (N_panels-1)

# Optical efficiency and power per ray
A_layout = panel_length * (panel_width*N_panels + spacing*(N_panels-1)) # Gross area including spacings
A_aperture = len(panel_positions) * panel_width * panel_length          # Mirrors' surface


# Main loop, where for each hour of the days in the datetime_list, the position of the sun is found, and the angle of tilt
# of the panels. The cumulative soiling is applied to the panels which then are put through pysoltrace to find the resulting
# efficiency of the plant. 
for i, date in enumerate(datetime_list):

    # Finding the position of the sun in the sky
    elevation_deg = solar.get_altitude(lat,lon,date) 
    azimuth_deg = solar.get_azimuth(lat,lon,date) 
    zenith_deg = 90 - elevation_deg  
    elevation_rad, azimuth_rad, zenith_rad = (np.deg2rad(x) for x in [elevation_deg, azimuth_deg, zenith_deg])
    incidence_angles = []

    # The absorbtion loss is put back to 0 at the start of each month, 
    # or after the desired cleaning frequency.
    if (date.day == 1 or date.day % cleaning_frequency == 0) and date.hour == 0:
        a_loss = np.zeros(len(panel_positions))
    
    if elevation_deg > 0:
        day_night = "day"
       
        # Calculate the DNI from the date and elevation, and the sun position vector (Sz, Sx, Sy). 
        # The first two elements of the sun position vector are then normalised. 
        solar_radiation = radiation.get_radiation_direct(date, elevation_deg) 
        sun_position = np.array([-np.sin(azimuth_rad)*np.sin(zenith_rad), np.cos(zenith_rad), np.cos(azimuth_rad)*np.sin(zenith_rad)])
        sn = sun_position[0:2]/np.linalg.norm(sun_position[0:2]) # sn is the NORMALISED sun position vector

        # Work needed here:
        # The angle of projection of the sun vector to the transverse plane
        theta_T = np.arctan(sn[0]/sn[1]) 
        theta_L = np.arctan(np.cos(azimuth_rad) * np.tan(zenith_rad))     
        r_aim = np.c_[-panel_positions,(receiver_height-mirrors_height)*np.ones_like(panel_positions)]                                       
        theta_aim = np.arctan(r_aim[:,0]/r_aim[:,1])
        tilt_deg = 180/np.pi*(theta_T+theta_aim)/2 
        tilt_rad = np.deg2rad(tilt_deg)
 
        sun_vector = np.transpose(np.array([np.cos(azimuth_rad)*np.cos(elevation_rad), 
                                            np.sin(azimuth_rad)*np.cos(elevation_rad), np.sin(elevation_rad)]))
        incidence_angle_rad = np.transpose(0.5*np.arccos(sun_vector.dot(towerToMirror_norm)))

        # Soiling Factor
        g_factor = 2/np.cos(incidence_angle_rad)          # Geometrical factor which takes into account blocking and shading effects.
        SR_hourly = 0.00042 * np.cos(tilt_rad)            # Constant Soiling Factor.
        SR_cumulated = SR_hourly + a_loss                 # Taking into account the buildup of dust already present.
        reflectivity_actual = reflectivity * (1-(SR_cumulated * g_factor))

        shift_x_axis = (stg1_height-mirrors_height)/np.tan(elevation_rad) # east-west axis [x-axis] --> negative shift (west is positive)
        shift_z_axis = (stg1_height-mirrors_height)*np.tan(theta_L)       # north-south [z-axis] ) --> positive shift (north is positive)

        break

print(tilt_deg)
