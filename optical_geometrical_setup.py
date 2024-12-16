# Optical and Geometrical Setup for PySolTrace

from pysoltrace import Point

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

def op_fictitious_surface(PT, slope_error, specularity_error):
    optics_fictitious = PT.add_optic("Fictitious Surface")
    optics_fictitious.front.reflectivity = 0.0 
    optics_fictitious.back.reflectivity = 0.0
    optics_fictitious.front.transmissivity = 1.0
    optics_fictitious.back.transmissivity = 1.0
    optics_fictitious.front.refraction_real = 1.0
    optics_fictitious.back.refraction_real = 1.0 
    optics_fictitious.front.slope_error = slope_error
    optics_fictitious.front.specularity_error = specularity_error
    optics_fictitious.back.slope_error = slope_error
    optics_fictitious.back.specularity_error = specularity_error
    return optics_fictitious

def op_cover_surface(PT, slope_error, specularity_error):
    optics_secondary = PT.add_optic("Cover")
    optics_secondary.front.reflectivity = 1.0 
    optics_secondary.front.transmissivity = 0.0
    optics_secondary.back.reflectivity = 0.0 
    optics_secondary.back.transmissivity = 0.0
    optics_secondary.front.slope_error = slope_error
    optics_secondary.front.specularity_error = specularity_error
    optics_secondary.back.slope_error = slope_error
    optics_secondary.back.specularity_error = specularity_error
    return optics_secondary

def op_heliostat_surface(PT, slope_error, specularity_error, reflectivity, p):
    optics_p = PT.add_optic(f"Heliostat_{p}") 
    optics_p.front.reflectivity = reflectivity  
    optics_p.front.transmissivity = 0.0
    optics_p.front.slope_error = slope_error
    optics_p.front.specularity_error = specularity_error
    optics_p.back.reflectivity = 0.0
    optics_p.back.transmissivity = 0.0
    return optics_p

def op_receiver_surface(PT):
    optics_receiver = PT.add_optic("Receiver")
    optics_receiver.front.reflectivity = 0.06
    optics_receiver.back.reflectivity = 0.06
    optics_receiver.front.transmissivity = 0.0
    optics_receiver.back.transmissivity = 0.0
    return optics_receiver

def op_secondaryReflector_surface(PT, slope_error, specularity_error):
    optics_secondary = PT.add_optic("Secondary Reflector")
    optics_secondary.front.reflectivity = 0.9
    optics_secondary.front.transmissivity = 0.0
    optics_secondary.back.reflectivity = 0.9
    optics_secondary.back.transmissivity = 0.0
    optics_secondary.front.slope_error = slope_error
    optics_secondary.front.specularity_error = specularity_error
    optics_secondary.back.slope_error = slope_error
    optics_secondary.back.specularity_error = specularity_error
    return optics_secondary

def trapezoidal_secondary_reflector(stg, optics, receiver_height, receiver_length):

    spos1 = [-0.15, 0, receiver_height-0.06]
    spos2 = [0.15, 0, receiver_height-0.06]
    spos3 = [0, 0, receiver_height+0.04]

    el = stg.add_element() 
    el.optic = optics
    el.position = Point(*spos1) 
    el.aim = Point(-300, 0, 275) 
    el.aperture_rectangle(0.27, receiver_length); 
    el.surface_flat(); 
    el.interaction = 2 

    el = stg.add_element() 
    el.optic = optics
    el.position = Point(*spos2) 
    el.aim = Point(300, 0, 275) 
    el.aperture_rectangle(0.27, receiver_length)
    el.surface_flat(); 
    el.interaction = 2 

    el = stg.add_element() 
    el.optic = optics
    el.position = Point(*spos3)  
    el.aim = Point(0.000001,0,0)
    el.aperture_rectangle(0.12, receiver_length)
    el.surface_flat()
    el.interaction = 2

    return el
