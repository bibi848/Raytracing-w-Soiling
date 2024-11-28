
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
start_date = datetime(2023, 1, 1,hour=0,minute=0,second=0,tzinfo=dt.timezone(timezoneOffset)) # Define the start date
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
       
        # Calculate the DNI from the date and elevation, and the sun position vector (Sz, Sy, Sx).
        # The sun_position vector is calculated in the aiming strategy document, along with the rest of the 
        # maths presented below.
        solar_radiation = radiation.get_radiation_direct(date, elevation_deg) 
        sun_position = np.array([np.sin(azimuth_rad)*np.sin(zenith_rad), np.cos(zenith_rad), np.cos(azimuth_rad)*np.sin(zenith_rad)])
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

        # Run Soltrace and add the sun
        PT = PySolTrace()
        sun = PT.add_sun()
        sun.position.x = 1000*sun_position[0] # Can maybe optimise this
        sun.position.y = 1000*sun_position[1]
        sun.position.z = 1000*sun_position[2]

        # OPTICAL PROPERTIES
        # Fictitious surface, Stage 1
        optics_fictitious = stg_surface(PT, "Fictitious Surface", 0.0, 0.0, 1.0, 1.0, slope_error, specularity_error)
        optics_fictitious.front.refraction_real = 1.0
        optics_fictitious.back.refraction_real = 1.0 

        # Cover, Stage 2
        optics_secondary = stg_surface(PT, "Cover", 1.0, 1.0, 0.0, 0.0, slope_error, specularity_error)

        # Mirrors, Stage 3
        optics_helios = []
        for p, x in enumerate(panel_positions):
            optics_p = stg_surface(PT, f"Heliostat_{p}", reflectivity_actual[p], 0.0, 0.0, 0.0, slope_error, specularity_error)
            optics_helios.append(optics_p)
        
        # Receiver, Stage 4
        optics_receiver = stg_surface(PT, "Tube", reflectivity_rec, reflectivity_rec, 0.0, 0.0)

        # Secondary reflector, Stage 4
        optics_secondary = stg_surface(PT, "Secondary", 0.9, reflectivity[0], 0.0, 0.0, slope_error, specularity_error)

        # GEOMETRY
        # Fictitious Surface, Stage 1
        st1 = PT.add_stage()
        st1.is_multihit = True  # True = flagged ; False = Unmarked 
        st1.is_virtual = False
        st1.name = 'Sdage 1 - Fictitious'
        st1.position = Point(0, 0, 0)
        el1 = st1.add_element()
        
        if azimuth_deg > 0 and azimuth_deg < 20:  
            if elevation_deg > 50:             
                el1.position = Point(-shift_x_axis/2.2 + 1e-06, stg1_height, shift_z_axis + 1e-06) 
            else:
                el1.position = Point(-0.8 + 1e-06, stg1_height, shift_z_axis + 1e-06) 
            el1.aperture_rectangle(stg1_length, stg1_width)
            el1.aim = Point(np.abs(shift_x_axis/2.2), 1000, shift_z_axis)
        
        if azimuth_deg > 20 and azimuth_deg < 90:  
            if elevation_deg >  30:
                el1.position = Point(-shift_x_axis/2.2 + 1e-06, stg1_height, shift_z_axis + 1e-06) 
            else: 
                el1.position = Point(-shift_x_axis/1.2 + 1e-06, stg1_height, shift_z_axis + 1e-06)
            el1.aperture_rectangle(stg1_length, stg1_width)
            el1.aim = Point(np.abs(shift_x_axis/2.2), 1000, shift_z_axis)
                    
        if azimuth_deg > 90 and azimuth_deg < 180:
            if elevation_deg > 50:
                el1.position = Point(-shift_x_axis/2.2 + 1e-06, stg1_height, shift_z_axis + 1e-06)
            else:
                el1.position = Point(-shift_x_axis + 1e-06, stg1_height, shift_z_axis + 1e-06)
            el1.aperture_rectangle(stg1_length, stg1_width)
            el1.aim = Point(np.abs(shift_x_axis/2.2), 1000, shift_z_axis)
                    
        if azimuth_deg > 180 and azimuth_deg < 270:
            if elevation_deg > 50:
                el1.position = Point(shift_x_axis/2.2 + 1e-06, stg1_height, shift_z_axis + 1e-06) 
            else:   
                el1.position = Point(shift_x_axis + 1e-06, stg1_height, shift_z_axis + 1e-06)
            el1.aperture_rectangle(stg1_length, stg1_width)
            el1.aim = Point(np.abs(shift_x_axis/2.2), 1000, shift_z_axis)
                    
        if azimuth_deg > 270 and azimuth_deg < 349: 
            if elevation_deg > 50:
                el1.position = Point(shift_x_axis/2.2 + 1e-06, stg1_height, shift_z_axis + 1e-06) 
            else:
                el1.position = Point(shift_x_axis + 1e-06, stg1_height, shift_z_axis + 1e-06)
            el1.aperture_rectangle(stg1_length, stg1_width)
            el1.aim = Point(-shift_x_axis/2.2, 1000, shift_z_axis)    

        if azimuth_deg > 349:
            if elevation_deg > 50:
                el1.position = Point(shift_x_axis/2.2 + 1e-06, stg1_height, shift_z_axis + 1e-06)
            else: 
                el1.position = Point(0.8 + 1e-06, stg1_height, shift_z_axis + 1e-06)
            el1.aperture_rectangle(stg1_length, stg1_width)
            el1.aim = Point(-shift_x_axis/2.2, 1000, shift_z_axis) 
            
        el1.surface_flat()
        el1.interaction = 1           # 2 for reflection, 1 for refraction 
        el1.optic = optics_fictitious # this shouldn't matter, but we need it? <- From Federico

        # Cover, Stage 2
        st2 = PT.add_stage()
        st2.is_multihit = True         
        st2.is_tracethrough = True
        st2.position = Point(0, 0, 0) 
        st2.name = 'Stage 2 - Cover'

        # 1st diagonal of trapezoid
        el2 = st2.add_element() 
        el2.optic = optics_secondary
        spos1 = [-0.15, receiver_height-0.06, 0] 
        el2.position = Point(*spos1) 
        el2.aim = Point(-300,275,0) 
        el2.aperture_rectangle(receiver_length, 0.27); 
        el2.surface_flat(); 
        el2.interaction = 2 

        # 2nd diagonal of trapezoid
        el2 = st2.add_element() 
        el2.optic = optics_secondary
        spos2 = [0.15,receiver_height-0.06,0] 
        el2.position = Point(*spos2) 
        el2.aim = Point(300,275,0) 
        el2.aperture_rectangle(receiver_length,0.27)
        el2.surface_flat(); 
        el2.interaction = 2 

        # Small base of trapezoid
        el2 = st2.add_element() 
        el2.optic = optics_secondary
        spos3 = [0,receiver_height+0.04,0]
        el2.position = Point(*spos3)  
        el2.aim = Point(0.000001,0,0)
        el2.aperture_rectangle(receiver_length,0.12)
        el2.surface_flat()
        el2.interaction = 2

        # Mirrors, Stage 3
        st3 = PT.add_stage()
        st3.is_multihit = True
        st3.position = Point(0, 0, 0) 
        st3.name = 'Stage 3 - Solar field'     
        
        # Separate evaluation of Xpos and Xaim because if they are the same value, SolTrace will fail the 
        # mirror's orientation. A small tolerance value has been added. 
        for i in range(len(panel_positions)): 
            el3 = st3.add_element()
            el3.optic = optics_helios[i]
            n = panel_normal(tilt_rad[i]) 
            hpos = np.array([panel_positions[i], mirrors_height, 0])
            aim = hpos + 1000*n 
            el3.position = Point(*hpos)          

            hpos = np.array([0.000001+panel_positions[i],0,0]) 
            aim = hpos + 1000*n                                
            el3.aim = Point(*aim)
            el3.aperture_rectangle(panel_length,panel_width)
            el3.surface_flat()
            el3.interaction = 2 

        # Secondary Reflector, Stage 4
        st4 = PT.add_stage()
        st4.name = 'stage 4 - rec + CPC'
        st4.position = Point(0, 0, 0)
        st4.is_multihit = True

        # Receiver, Stage 4
        el4 = st4.add_element()
        el4.optic = optics_receiver
        rpos = [0,receiver_height,0]
        el4.position = Point(*rpos)
        el4.aim = Point(0,0,0)
        el4.aperture_singleax_curve(0,0,receiver_length)
        el4.surface_cylindrical(receiver_diameter/2.0)
        el4.interaction = 2

        # Secondary geometry (trapezoidal geometry)
        el4 = st4.add_element() 
        el4.optic = optics_secondary
        spos1 = [-0.15,receiver_height-0.07,0] 
        el4.position = Point(*spos1) 
        el4.aim = Point(-300,275,0) 
        el4.aperture_rectangle(receiver_length,0.27); 
        el4.surface_flat(); 
        el4.interaction = 2 

        # 2nd trapezoidal diagonal
        el4 = st4.add_element() 
        el4.optic = optics_secondary
        spos2 = [0.15,receiver_height-0.07,0] 
        el4.position = Point(*spos2) 
        el4.aim = Point(300,275,0) 
        el4.aperture_rectangle(receiver_length,0.27);
        el4.surface_flat(); 
        el4.interaction = 2 

        # Small base of trapezoid
        el4 = st4.add_element() 
        el4.optic = optics_secondary
        spos3 = [0,receiver_height+0.03,0]
        el4.position = Point(*spos3)  
        el4.aim = Point(0.000001,0,0)
        el4.aperture_rectangle(receiver_length,0.12)
        el4.surface_flat()
        el4.interaction = 2
            
        PT.num_ray_hits = number_hits 
        PT.max_rays_traced = PT.num_ray_hits*100 
        PT.is_sunshape = sunshape_flag          # True
        PT.is_surface_errors = sfcerr_flag      # True
        PT.dni= 1000

        # RAYTRACING   
        PT.run(-1, False)
        df = PT.raydata # Dataframe that collects all the data related to the simulation

        # Add one more column at the actual dataframe in order to simplify the further analysis
        counts = df['number'].value_counts() 
        counts_dict = counts.to_dict()                    # Convert the counts series to a dictionary
        df['occurrences'] = df['number'].map(counts_dict) # Create a new column 'occurrences' based on the counts
                
        # Cover analysis 
        cover_df = df[(df['stage'] == 2) & (df['element'] != 0)]['number'].unique().shape[0] # Absorbed or reflected by the cover
        cover_miss = df[(df['stage']==2) & (df['element'] == 0)]['number'].unique().shape[0] 
        rays_stage2 = cover_df + cover_miss

        # Mirrors analysis
        mirrors_abs = df[(df['stage']==3) & (df['element'] < 0)]['number'].unique().shape[0]     # Reflectance loss mirrors
        mirrors_refl = df[(df['stage'] == 3) & (df['element'] > 0)]['number'].unique().shape[0]  # Rays may or may not hit the stage 3
        mirrors_hits = df[(df['stage']==3) & (df['element'] != 0)]['number'].unique().shape[0]   # Can be reflected, absorbed, or miss,
                                                                                                 # going in the gap between mirrors
        rays_gaps = cover_miss - mirrors_hits                # Rays in the space between mirrors
        rays_stage3 = mirrors_abs + mirrors_refl + rays_gaps  

        # Receiver analysis
        receiver_abs = df[(df['stage'] == 4) & (df['element'] == -1)].shape[0] # number of rays hit the target and have been reflected by the mirrors
        receiver_refl =  df[(df['stage'] == 4) & (df['element'] == 1)].shape[0]
        receiver_tot = receiver_abs + receiver_refl

        # cosine losses
        A_eff_i = panel_length * (panel_width * np.cos(tilt_rad))  # Effective surface for each panel due to slope angle
        A_field = np.sum(A_eff_i)  
        A_tilted = A_field  
        avg_cos_beta = np.mean(np.cos(tilt_rad))

        # Efficiency
        ppr_corrected = total_width * panel_length * np.cos(theta_T) * PT.dni / (number_hits-rays_gaps)

        if (mirrors_hits) != 0:
            eta_opt_corrected = receiver_abs * ppr_corrected / (PT.dni * A_aperture)
        else: eta_opt_corrected = 0
        
        E_sun = A_aperture * solar_radiation       # overall amount of energy that hits the total aperture area
        thermal_energy = eta_opt_corrected * E_sun # thermal energy at the receiver
        E_sun_real = E_sun * eta_opt_corrected

    else:
        day_night = "night"
        eta_opt_corrected = 0
        SR_cumulated = a_loss #during nigth hours there is no increase in dust (it keeps constant)
        reflectivity_actual = reflectivity * (1-SR_cumulated)
        # theta_i = np.zeros_like(panel_positions) 
        # to model the vertical position of panels during night hours
        # NO DUST IS DEPOSITED over the mirrors' surface    
        incidence_angle_rad = np.full_like(panel_positions, np.nan)  
        g_factor = np.zeros_like(panel_positions) 
        tilt_deg = 90 * np.ones_like(panel_positions)
        solar_radiation = 0
        E_sun = 0
        E_sun_real = 0
        thermal_energy = 0
    
    a_loss = SR_cumulated # to update the value of a_loss for the next iteration
    eta_opt_zero = 0.686
    IAM = eta_opt_corrected/eta_opt_zero   
    month = date.month
    day = date.day
    hour = date.hour
    rho_avg_field = np.mean(reflectivity_actual) # 1 value for the whole plant 
    
    lookup_table_soiled.loc[lookup_table_loc] = [date, azimuth_deg, elevation_deg, eta_opt_corrected, IAM, 
                                             day_night, month, day, hour, rho_avg_field, solar_radiation, 
                                             E_sun, E_sun_real, thermal_energy] + list(np.rad2deg(incidence_angle_rad)) + list(tilt_deg) + list(reflectivity_actual) + list(a_loss) + list(g_factor)
    lookup_table_loc += 1

print(lookup_table_soiled)
timer_end = time.time()

print(timer_end - timer_start)

lookup_table_soiled.to_csv(file_path, index=False)
