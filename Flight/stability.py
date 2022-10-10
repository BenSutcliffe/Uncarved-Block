import numpy as np
from scipy.integrate import quad

def CNalphaN_subs(num_fin, fin_span, Area_ref, Area_fin, Beta, mid_chord_gamma):
  """
    Computes the normal force coefficient for all the fins in the subsonic regime, M<0.8.
    Barrowman analysis

    Parameters
    ----------
    num_fin : float
      number of fins on the rocket
    fin_span : float
      fin span height, distance from root to tip
    Area_ref: float
      cross sectional reference area for the rocket body
    Area_fin: float
      fin planar side area
    Beta: float
      root(Machno^2 - 1), inverse of Prandtl number
    mid_chord_gamma: float
      angle between the two midchord positions of the fins
      
    Returns
    -------
    CNalphaN: float
      value of the normal force coefficient of the fins in subsonic regime
  """
  CNalphaN = (num_fin/2) * (2 * np.pi) * (((fin_span)**2) /Area_ref)/(1 + np.sqrt(1 + ((Beta * (fin_span)**2)/(Area_fin * np.cos(mid_chord_gamma)))**2))
  return CNalphaN

def CNalphaN_super(num_fin, Area_ref, Area_fin, Beta, angle_attack):
  """
    Computes the normal force coefficient for all the fins in the supersonic regime M>1.2 based
    on Busemann theory.

    Parameters
    ----------
    num_fin : float
      number of fins on the rocket
    Area_ref: float
      cross sectional reference area for the rocket body
    Area_fin: float
      fin planar side area
    Beta: float
      root(Machno^2 - 1), inverse of Prandtl number
    angle_attack: float
      angle of attack of the vehicle
      
    Returns
    -------
    CNalphaN: float
      value of the normal force coefficient of the fins in supersonic regime
  """
  #Computes the Mach value corresponding with beta
  M_value = np.sqrt(np.abs(1-(Beta)**2))

  #Finds the values for the three velocity dependent coefficients
  K_1 = 2/Beta 
  K_2 = (((2.4*M_value**4)-(4*Beta**2))/(4*Beta**4))
  K_3 = (((2.4*M_value**8)-(10.88*M_value**6)+(24*M_value**4)+(8))/(6*Beta**7))

  #Finds the values for the force coefficient using the above coefficients and angle of attack
  CNalphaN = (num_fin/2) * (Area_fin/Area_ref) * ((K_1) + (K_2 * angle_attack * np.pi/180) + (K_3 * (angle_attack * np.pi/180)**2)) 
  return CNalphaN

def CNalphaN_trans(num_fin, fin_span, Area_ref, Area_fin, Mach_no, mid_chord_gamma, angle_attack):
  """
    Linearly interpolates normal force coefficient for the fins for 0.8<M<1.2

    Parameters
    ----------
    num_fin : float
      number of fins on the rocket
    fin_span : float
      fin span height, distance from root to tip
    Area_ref: float
      cross sectional reference area for the rocket body
    Area_fin: float
      fin planar side area
    Mach_no: float
      Mach number of the vehicle
    mid_chord_gamma: float
      angle between the two midchord positions of the fins
      
    Returns
    -------
    CNalphaN: float
      value of the normal force coefficient of the fins in transsonic regime
  """
  beta_sub = 0.6 #Value of inverse of Prandtl number for M=0.8
  beta_super = 0.6633
  CNalphaNpoint8 = CNalphaN_subs(num_fin, fin_span, Area_ref, Area_fin, beta_sub, mid_chord_gamma)
  CNalphaNonetwo = CNalphaN_super(num_fin, Area_ref, Area_fin, beta_super, angle_attack)
  CNalphaN = CNalphaNpoint8 + ((CNalphaNonetwo-CNalphaNpoint8)*((Mach_no-0.8)/(0.4)))
  return CNalphaN

def C_N(Area_plan, Area_ref, angle_attack, K=1.1):
  """
    Returns the normal force coefficient for the body tube

    Parameters
    ----------
    Area_plan : float
      side area of the body
    Area_ref: float
      cross sectional reference area for the rocket body
    Area_fin: float
      fin planar side area
    angle_attack: float
      angle of attack of the vehicle
    K:
      constant scalar term
      
    Returns
    -------
    CNalphaN: float
      value of the normal force coefficient of the body tube
  """
  C_N = K * Area_plan/Area_ref * (np.sin(np.pi*angle_attack/180))**2
  return C_N

def X_N(nose_height, body_height):
  """
    Returns the centre of pressure of the body tube relative to the top of the nosecone

    Parameters
    ----------
    nose_height : float
      height of the nosecone
    body_height: float
      height of the body tube

      
    Returns
    -------
    CP_bodytube: float
      vertical position of the centre of pressure of the body tube relative to the top of the nosecone
  """
  CP_bodytube = body_height/2 + nose_height #Adds half the height to the nosecone length
  return CP_bodytube

def X_sf(c, s, A_fin, beta):
  """
    Returns the centre of pressure of the fins relative to their leading point for the supersonic regime

    Parameters
    ----------
    nose_height : float
      height of the nosecone
    body_height: float
      height of the body tube

      
    Returns
    -------
    CP_bodytube: float
      vertical position of the centre of pressure of the body tube relative to the top of the nosecone
  """
  #Method to find the new COP for the fins in the supersonic regime
  A_r = 2*(s*0.01)**2 / A_fin #Calculates the fin aspect ratio
  X_f = c * (((A_r * beta) - 0.67)/((2*A_r*beta) - 1)) #Finds the COP of the fin, as in OpenRocket
  return X_f

def c_g(x, Fins):
  #Equation for defining chord with fin height
    c_sq = (Fins.Chord_root + ((Fins.Chord_tip - Fins.Chord_root)*x)/Fins.fin_span)**2
    return c_sq

def MAC(Fins, c_f):
  #Integration to find MAC length
  result, err = quad(c_f,0,Fins.fin_span,args=(Fins,))
  Mac = result/(Fins.area() * 10000)
  return Mac

def c_LE(x, Fins):
  #Product of chord equation and LE equation
    c_sq = (Fins.Chord_root + ((Fins.Chord_tip - Fins.Chord_root)*x)/Fins.fin_span)*((Fins.Chord_root - Fins.Chord_tip)*x/Fins.fin_span)
    return c_sq

def MAC_x(Fins, c_func):
  #Integration to find MAC x position from edge
  result, err = quad(c_func,0,Fins.fin_span,args=(Fins,))
  Mac = result/(Fins.area() * 10000)
  return Mac

def pressure_position_transonic(c, s, A_fin, Mach):
  f_1 = X_sf(c, s, A_fin, 1.732)/c
  dB = 1e-2
  beta_1 = np.sqrt(np.abs((2+dB)**2 - 1))
  beta_2 = np.sqrt(np.abs((2-dB)**2 - 1))
  f_2 = (X_sf(c, s, A_fin, (beta_1))/c  - X_sf(c, s, A_fin, (beta_2))/c)/(2*dB)
  a = np.array([[1/32, 1/16, 1/8, 1/4, 1/2, 1], [5/16, 1/2, 3/4, 1, 1, 0], [32, 16, 8, 4, 2, 1], [80, 32, 12, 4, 1, 0], [160, 48, 12, 2, 0, 0], [240, 48, 6, 0, 0, 0]])
  b = np.array([0.25, 0, f_1, f_2,  0, 0])
  x = np.linalg.solve(a, b)
  value = ((x[0] * (Mach**5)) + (x[1] * (Mach**4)) + (x[2] * (Mach**3)) + (x[3] * (Mach**2)) + (x[4] * (Mach)) + (x[5]))*c
  return value

def F_fin_N(C_N_, rho, A, v, alpha):
  F_fin_val = 0.5 * C_N_ * rho * A * (v**2)
  return F_fin_val

def C_N_force(C_N_a, delta):
  C_m_val = C_N_a*delta
  return C_m_val

def M_fin_norm(F_N, X):
  M_fin_n = F_N * (X)
  return M_fin_n