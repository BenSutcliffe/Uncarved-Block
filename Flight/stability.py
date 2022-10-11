import numpy as np
from scipy.integrate import quad

G_alu = 24E9
P_atm = 1E5

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

def CNalphaN_trans(n, s, A_ref, A_fin, Mach, gamma, alpha):
  CNalphaNpoint8 = CNalphaN_subs(n, s, A_ref, A_fin, 0.6, gamma)
  CNalphaNonetwo = CNalphaN_super(n, A_ref, A_fin, 0.663, alpha)
  CNalphaN = CNalphaNpoint8 + ((CNalphaNonetwo-CNalphaNpoint8)*((Mach-0.8)/(0.4)))
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

def flutter(G, h, c_r, c_t, M, P):
    s = 0.5*(c_r + c_t)*h
    a = ((h**2)/s)
    lam = c_t/c_r
    t = (((M**2) * (P * (lam + 1)*(a**3))/(1.5*G*(a+2)))**(1/3))*c_r
    return t

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