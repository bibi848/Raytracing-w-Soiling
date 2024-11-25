
# Commercial Trapezoidal 4 Stages Lookup Table Soiled, 14 Day Cleaning Rotation

# System Setup
# Coordinate System, XYZ: X = West, Y = Perpendicular to the ground, Z = North.
# The 1st mirror is the East-most mirror, and the 11th mirror is the West-most mirror.

# FUNCTIONS

import numpy as np
import pysolar.solar as solar
import pysolar.radiation as radiation
from pysoltrace import PySolTrace, Point
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import datetime as dt
from IPython.display import display, Latex
import st_processing_functions
import ipywidgets as widgets
from IPython.display import display
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.optimize import bisect
from datetime import datetime, timedelta
import pytz
import calendar
from pvlib import irradiance
import pvlib


# Location: Woomera
lat, lon = lat,lon = -31.2,136.816667
timezoneOffset = dt.timedelta(hours = 9.5)

# Time Duration: Measuring across a whole year
# It starts on the 1st Jan 2024, and ends on the 1st Jan 2025, with a date for every hour of the year.
start_date = datetime(2024, 1, 1,hour=0,minute=0,second=0,tzinfo=dt.timezone(timezoneOffset)) # Define the start date
end_date = start_date + timedelta(days=365) # Calculate the end date 
datetime_list = []                          # Create a list of datetime objects
current_date = start_date                   # Loop to generate datetime objects from start_date to end_date
while current_date <= end_date:
    datetime_list.append(current_date)
    current_date += timedelta(hours=1)

# Simulation Parameters
number_hits = 1e3  # 1 million rays 
sunshape_flag = False
sfcerr_flag = False

# Geometry of Panels
total_width = 7.5    #m
spacing = 0.2        #m
panel_width = 0.5    #m
panel_length = 12.18 #m
mirrors_height = 0.5 #m
panel_positions = np.arange(-3.5, 3.75, panel_width + spacing)
N_panel = len(panel_positions) # Number of panels, 11














