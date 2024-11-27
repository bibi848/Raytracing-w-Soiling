
# Optical and Geometrical Setup for PySolTrace

import numpy as np
import pysolar.solar as solar
import pysolar.radiation as radiation
from pysoltrace import PySolTrace, Point
import datetime as dt
import pandas as pd


# Optical
def stg_surface(PT, name, f_reflectivity, b_reflectivity, f_transmissivity, b_transmissivity, slope_error = 0, specularity_error = 0):
    optics = PT.add_optic(name)
    optics.front.reflectivity = f_reflectivity
    optics.back.reflectivity = b_reflectivity
    optics.front.transmissivity = f_transmissivity
    optics.back.transmissivity = b_transmissivity
    optics.front.slope_error = slope_error
    optics.front.specularity_error = specularity_error
    optics.back.slope_error = slope_error
    optics.back.specularity_error = specularity_error
    return optics