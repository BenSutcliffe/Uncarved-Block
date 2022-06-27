import numpy as np

def flutter(G, h, c_r, c_t, t, P):
    s = 0.5*(c_r + c_t)*h
    a = (h**2)/s
    lam = c_t/c_r
    flutter_square = (1.5*G*(a+2)*((t/c_r)**3))/(P * (lam + 1)*(a**3))
    flutter_mach = np.sqrt(flutter_square)
    mass = 0.5*(c_r+c_t)*h*t*2710*4
    return (flutter_mach, mass)

g_1 = 24E9
p_1 = 1E5
c_r1 = 0.8
c_t1 = 0.6
h1 = 0.33

t1 = 0.01

print(flutter(g_1, h1, c_r1, c_t1, t1, p_1))