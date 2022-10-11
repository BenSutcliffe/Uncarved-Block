"""
Stability Function - Calculates the stability of the vehicle
"""
#Import all the required functions and classes from the other modules
import numpy as np
from aero_coefficients import MAC_length, chord_sq, MAC_x_distance, c_LE
from Global_vars import *
from aero_coefficients import Xf_supersonic, X_N, C_N_body, CNalphaN_subs, CNalphaN_super, Xf_transonic
import Classes

#Imports Body and nosecone objects from classes
A_body = Classes.Bodyone
vehicle_nosecone = Classes.Cone

def Griffin_stability(Mach_number, vehicle_length, nosecone, vehicle_fins, vehicle_body, angle_attack, number_fins, COM):
    """
    Returns the stability of the vehicle given a specific Mach number, finset, body tube and nosecone. Allows user to
    see how stability changes with Mach number.

    Parameters
    ----------
    Mach_number : float
      mach number at which the stability of the vehicle is desired
    vehicle_length : float
      total vehicle length
    ambient_pressure : float
        external pressure at maximum Q, static
    nosecone : Class
        nosecone object for the vehicle
    vehicle_fins : Class
        Fin object for the vehicle
    vehicle_body : Class
        Body object for the vehicle
    angle_attack : float
        angle of attack of the vehicle at the test point
    number_fins : float
        number of fins on the vehicle
    COM : float
        centre of mass of the vehicle
      
    Returns
    -------
    stability_Cal : float
      stability of the vehicle in calibres
    """
    
    beta = np.sqrt(np.abs(1-(Mach_number)**2)) #Calculates the inverse Prandtl number
    fin_pos_Panthera = vehicle_length - vehicle_fins.Chord_root #finds the top of the fins relative to the nosecone tip

    #Finds the fin centre of pressure & normal force coefficent based on Mach number
    if Mach_number <= 0.8: 
        #Totally subsonic Regime
        CNangle_attack = CNalphaN_subs(number_fins, vehicle_fins.fin_span, vehicle_body.Arearef(), vehicle_fins.area(), beta, vehicle_fins.fin_gamma())
        CP_x_fins = (vehicle_fins.X_f() + fin_pos_Panthera) # Finds the centre of pressure for finset relative to the top of the rocket

    elif Mach_number < 1.2:
        #fin CP position shifts for subsonic -> transonic
        CNangle_attack = CNalphaN_subs(number_fins, vehicle_fins.fin_span, vehicle_body.Arearef(), vehicle_fins.area(), beta, vehicle_fins.fin_gamma())
        CP_x_fins = Xf_transonic(MAC_length(vehicle_fins, chord_sq), vehicle_fins.fin_span, vehicle_fins.area(), Mach_number) + fin_pos_Panthera + MAC_x_distance(vehicle_fins, c_LE)
    
    elif Mach_number <= 2:
        #normal force coefficient changes to supersonic formula, CP position continues to move transonic -> supersonic
        CNangle_attack = CNalphaN_super(number_fins, vehicle_body.Arearef(), vehicle_fins.area(), beta, angle_attack)
        CP_x_fins = Xf_transonic(MAC_length(vehicle_fins, chord_sq), vehicle_fins.fin_span, vehicle_fins.area(), Mach_number) + fin_pos_Panthera + MAC_x_distance(vehicle_fins, c_LE)

    else:
        #Pure supersonic, CP position shifts under supersonic formula
        CNangle_attack = CNalphaN_super(number_fins, vehicle_body.Arearef(), vehicle_fins.area(), beta, angle_attack)
        CP_x_fins = (Xf_supersonic(MAC_length(vehicle_fins, chord_sq), vehicle_fins.fin_span, vehicle_fins.area(), beta) + fin_pos_Panthera + MAC_x_distance(vehicle_fins, c_LE))

    CP_x_nosecone = nosecone.X_f() # Finds the centre of pressure for the nosecone relative to the nosecone tip
    CP_x_body = X_N(nosecone.nose_height, (vehicle_body.body_height)) # Finds the centre of pressure for body relative to the nosecone tip

    Normal_coeff_nosecone = nosecone.CombinedC #Normal force coefficent nosecone
    Normal_coeff_fins = CNangle_attack * vehicle_fins.K() #Finds the normal force coefficent for a finset with interference
    Normal_coeff_body = C_N_body(vehicle_body.areat(), vehicle_body.Arearef(), angle_attack) #Finds the normal force coefficent for the body

    #Finds the overall centre of pressure of the vehicle
    overall_COP = ((CP_x_nosecone * Normal_coeff_nosecone) + (CP_x_body * Normal_coeff_body)+ (CP_x_fins * Normal_coeff_fins))/(Normal_coeff_nosecone + Normal_coeff_body + Normal_coeff_fins)

    stability_Cal = (overall_COP - COM)/vehicle_body.body_diameter #stability of the vehicle in calibres
    return stability_Cal