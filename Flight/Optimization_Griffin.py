"""
Optimization Solver - Run this file to generate the fin
"""
#Import all the required functions and classes from the other modules
import numpy as np
import Classes
from Classes import Fins, Body
from flutter import flutter
from Global_vars import *
from aero_coefficients import CNalphaN_super
from forces import C_N_force, F_fin_N
from stability import Griffin_stability
from scipy.optimize import differential_evolution

#Imports Body and nosecone objects from classes
test_body = Classes.Bodyone
vehicle_nosecone = Classes.Cone
fin_thickness = 0 #initialise thickness
fin_force = 0 #initialise fin force

def objective(v):
    """
    Function to which differential optimization is applied. In this case, the function is used to optimise the finset 
    for minimum mass.

    Parameters
    ----------
    v : list
      list containing root chord length, tip chord length & fin span in m in that order
      
    Returns
    -------
    output : float
      returns mass or very large number
    """
    global fin_thickness #imports global variables
    global fin_force

    chord_root, chord_tip, fin_span = v #Assigns componenets of v to relevant lengths
    root_ratio = chord_tip/chord_root #computes fin root ratio

    Aspect_ratio = (fin_span**2)/(0.5*(chord_root + chord_tip)*fin_span) #computes the fin aspect ratio
    test_fins = Fins(chord_root, fin_span, chord_tip, (chord_root - chord_tip), Body_dia) #Initialises a fin object for the tested configuration

    fin_thickness = flutter((max_q_velo*1.2), speed_sound, max_q_staticp, p_atm, chord_root, Aspect_ratio, G_alu_psi, root_ratio) #calculates required thickness from flutter

    mass_fins = fin_thickness * test_fins.area() * density_alu * N_fins #computes fin mass given geometry and thickness
    
    Mach_numbers = np.linspace(0.3, 5.5, 30) #mach numbers to check stability at
    for velocity in Mach_numbers:
        stability_calibre = Griffin_stability(velocity, total_length, vehicle_nosecone, test_fins, test_body, angle_attack, N_fins, CoM) #calculates stability for this config
        if stability_calibre < desired_stability or chord_tip > (chord_root-0.15): #if it is too unstable or poor geometry
            output = 10e9 #returns an infeasible value - won't appear as mass optimal
            return output
        else:
            pass
    output = mass_fins #if sufficiently stable, outputs mass, which is optimised

    Normal_coeff = CNalphaN_super(N_fins, test_body.Arearef(), test_fins.area(), max_q_beta, angle_attack_force) #computes normal coeff at Max Q
    Normal_coeff_alpha = C_N_force(Normal_coeff, (angle_attack_force* np.pi/180)) #computes normal coeff with angle of attack
    fin_force = F_fin_N(Normal_coeff_alpha, max_q_rho, test_body.Arearef(), max_q_velo) #computes peak fin force at max Q

    return output

chord_limits = [0.01, 1.5] #sets chord length limits for optimizer
span_limits = [0.01, 0.7] #sets chord length limits for optimizer

#Computes the differential optimization
bounds = [chord_limits, chord_limits, span_limits]
result = differential_evolution(objective, bounds)
print('Status : %s' % result['message'])
print('Total Evaluations: %d' % result['nfev'])

# evaluate solution
solution = result['x']
evaluation = objective(solution)

#prints the optimized values
print('')
print('Root Chord length /m: = {}'.format(solution[0]))
print('Tip Chord length /m: = {}'.format(solution[1]))
print('Fin Span /m: = {}'.format(solution[2]))
print('Total mass /kg: = {}'.format(evaluation))
print('')
print('Fin thickness /mm: = {}'.format(fin_thickness*1000))
print('Max Fin Force /kN: = {}'.format(fin_force/1000))
