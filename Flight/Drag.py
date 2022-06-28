import numpy as np

def Reynolds(length, velocity, viscosity = 1.5e-5):
    #Calculates Reynolds number
    R = (length * velocity)/viscosity
    return R

def R_crit(Roughness, Length):
    #Finds critical Reynolds number for turbulence
    R_crit = 51 * (Roughness/Length)**(-1.039)
    return R_crit

def friction_Re(Reynold, R_critical, Roughness, totatlen):
    #Finds friction coefficent based on Reynolds number
    Coeff_friction = 0
    if Reynold <= 1e4:
        Coeff_friction = 1.48e-2
    elif Reynold > 1e4 and Reynold < R_critical:
        Coeff_friction = ((1.5*np.log(Reynold))-5.6)**(-2)
    elif Reynold >= R_critical:
        Coeff_friction = 0.032 * ((Roughness/totatlen)**0.2)
    return Coeff_friction


def friction_vel(Coeff_friction, Machno, A_ref, fineness, A_body_wet, fin_thickness, A_fin_wet, MAC_base):
    #Computes overall friction coefficent based on velocity
    Skin_friction_coeff = 0
    if Machno < 0.8:
        Skin_friction_coeff = Coeff_friction * (1 - ((0.1) * (Machno ** 2)))
    elif Machno > 0.8:
        Supersonic_skin_friction_1 = Coeff_friction * (1 + ((0.15) * (Machno ** 2)))**(-0.58)
        Supersonic_skin_friction_2 = Coeff_friction * (1 + ((0.18) * (Machno ** 2)))**(-1)
        if Supersonic_skin_friction_1 > Supersonic_skin_friction_2:
            Skin_friction_coeff = Supersonic_skin_friction_1
        else:
            Skin_friction_coeff = Supersonic_skin_friction_2
    Coeff_friction_total = (Skin_friction_coeff/A_ref) * (((1 + (1/(2 * fineness))) * (A_body_wet)) + ((1 + (2 * fin_thickness/MAC_base)) * (8*A_fin_wet)))
    return Coeff_friction_total


def Body_pressure_drag(Machno, maximumdiameter, noseconelen):
    #Calculates body pressure drag
    eta = np.arctan(maximumdiameter/(2*noseconelen))
    Coeff_body_pressure = 0
    Coeff_body_pressure_subs = 0.8*(np.sin(eta))**2
    Coeff_body_pressure_M1 = np.sin(eta)
    Coeff_body_pressure_super = (2.1 * ((np.sin(eta))**2)) + ((np.sin(eta))/(2 * np.sqrt((Machno**2) - 1)))
    if Machno < 0.8:
        Coeff_body_pressure = Coeff_body_pressure_subs
    elif Machno >= 0.8 and Machno < 1:
        Coeff_body_pressure = Coeff_body_pressure_subs + ((Coeff_body_pressure_M1-Coeff_body_pressure_subs)*(Machno-0.8)/0.2)
    elif Machno >= 1 and Machno <= 1.3:
        Coeff_body_pressure = Coeff_body_pressure_M1 + ((Coeff_body_pressure_super-Coeff_body_pressure_M1)*(Machno-1)/0.3)
    elif Machno >= 1.3:
        Coeff_body_pressure = Coeff_body_pressure_super
    return Coeff_body_pressure

def fin_pressure_drag(Machno, Leadingangle_base):
    #Calculates fin pressure drag
    Coeff_fin_pressure_LE = 0  #Assuming a rounded leading edge
    if Machno < 0.9:
        Coeff_fin_pressure_LE = ((1 - (Machno**2))**(-0.417)) - 1
    elif Machno >= 0.9 and Machno <1:
        Coeff_fin_pressure_LE = 1 - (1.785 * (Machno-0.9))
    elif Machno >= 1:
        Coeff_fin_pressure_LE = 1.214 - (0.502/(Machno**2)) + (0.1095/(Machno**4))

    Coeff_fin_pressure_TE = 0
    if Machno < 1:
        Coeff_fin_pressure_TE = 0.5*(0.12 + (0.13 * (Machno**2)))
    elif Machno >= 1:
        Coeff_fin_pressure_TE = 0.5*(0.25/Machno)
    
    Coeff_fin_pressure_base = Coeff_fin_pressure_LE * (np.cos(Leadingangle_base * np.pi/180))**2
    Coeff_base_drag = Coeff_fin_pressure_TE
    return(Coeff_fin_pressure_base, Coeff_base_drag)
