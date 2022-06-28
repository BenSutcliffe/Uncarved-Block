import numpy as np
from scipy.integrate import quad

def CNalphaN_subs(n, s, A_ref, A_fin, Beta, gamma):
  #This finds the normal force coefficient for a given fin
  CNalphaN = (n/2) * (2 * np.pi) * (((s*0.01)**2) /A_ref)/(1 + np.sqrt(1 + ((Beta * (s*0.01)**2)/(A_fin * np.cos(gamma)))**2)) #From the open rocket documentation
  return CNalphaN

def CNalphaN_super(n, A_ref, A_fin, Beta, alpha):
  #New method to find the normal force coefficient for the non-subsonic regime
  M_value = np.sqrt(np.abs(1-(Beta)**2))

  #Finds the values for the three velocity dependent coefficients
  K_1 = 2/Beta 
  K_2 = (((2.4*M_value**4)-(4*Beta**2))/(4*Beta**4))
  K_3 = (((2.4*M_value**8)-(11.88*M_value**6)+(24*M_value**4)+(8))/(6*Beta**7))

  #Finds the values for the force coefficient using the above coefficients and angle of attack
  CNalphaN = (n/2) * (A_fin/A_ref) * ((K_1) + (K_2 * alpha * np.pi/180) + (K_3 * (alpha * np.pi/180)**2)) 
  return CNalphaN

def C_N(A_plan, A_ref, alpha, K=1.1):
  #Normal force coefficient of the rocket body
  C_N = K * A_plan/A_ref * (np.sin(np.pi*alpha/180))**2
  return C_N

def X_N(nose, h):
  #Finds centre of pressure for the cross sectional area
  X_N = h/2 + nose #Adds half the height to the nosecone length
  return X_N

def X_sf(c, s, A_fin, beta):
  #Method to find the new COP for the fins in the supersonic regime
  A_r = 2*(s*0.01)**2 / A_fin #Calculates the fin aspect ratio

  X_f = c * (((A_r * beta) - 0.67)/((2*A_r*beta) - 1)) #Finds the COP of the fin, as in OpenRocket
  return X_f


def X(bottom, planform, nose, alpha, CNalpha):
  #Function to find the overall centre of pressure for the rocket
  
  X_1 = (bottom.X_f() + 224) # Finds the centre of pressure for finset relative to the top of the rocket
  X_3 = X_N(planform.n_1, planform.h) # Finds the centre of pressure for body relative to the top of the rocket
  X_4 = nose.X_f() # Finds the centre of pressure for the nosecone relative to the top of the rocket

  C_1 = CNalpha * bottom.K() #Finds the normal force coefficent for a finset
  C_3 = C_N(planform.areat(), planform.Arearef(), alpha) #Finds the normal force coefficent for the body
  C_4 = nose.C[alpha] #Finds the normal force coefficent for the nosecone

  #Finds the overall centre of pressure as a weighted mean of component COP and force coefficients, as from OpenRocket Documentation
  X = ((X_1 * C_1) + (X_3 * C_3) + (X_4 * C_4))/(C_1 + C_3 + C_4)
  return X

def c_g(x, Fins):
  #Equation for defining chord with fin height
    c_sq = (Fins.C_r + ((Fins.C_t - Fins.C_r)*x)/Fins.s)**2
    return c_sq

def MAC(Fins, c_f):
  #Integration to find MAC length
  result, err = quad(c_f,0,Fins.s,args=(Fins,))
  Mac = result/(Fins.area() * 10000)
  return Mac

def c_LE(x, Fins):
  #Product of chord equation and LE equation
    c_sq = (Fins.C_r + ((Fins.C_t - Fins.C_r)*x)/Fins.s)*((Fins.C_r - Fins.C_t)*x/Fins.s)
    return c_sq

def MAC_x(Fins, c_func):
  #Integration to find MAC x position from edge
  result, err = quad(c_func,0,Fins.s,args=(Fins,))
  Mac = result/(Fins.area() * 10000)
  return Mac