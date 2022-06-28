import numpy as np

def flutter(G, h, c_r, c_t, M, P):
    s = 0.5*(c_r + c_t)*h
    a = ((h**2)/s)
    lam = c_t/c_r
    t = (((M**2) * (P * (lam + 1)*(a**3))/(1.5*G*(a+2)))**(1/3))*c_r
    return t

G_alu = 24E9
P_atm = 1E5