import numpy as np
from scipy.integrate import quad

def flutter(v, a, p, p0, c, A, G_el, lam):
  result = (((v/a)**2)*39.3*(p/p0)*((lam+1)/2)*(A**3))/(G_el*(A+2))
  t = c * ((result)**(1./3))
  return t