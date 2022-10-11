"""
Body flight profile critical values
"""

"""
Taken from Henry Free's preliminary design document - Flight parameters
"""
start_mass = 632 #Wet mass, kg
end_mass = 208 #Recovery Mass, kg
v_burnout = 1588.34 #Burnout velocity, m/s
max_q_velo = 841.079701 #Velocity at Max-Q, m/s
max_q_rho = 0.403219 #Density at Max-Q, kg/m3
max_q_staticp = 25687.78 #Static pressure at Max-Q, Pa
p_atm = 94321.68

"""
Taken from Henry Free's preliminary design document - Geometry parameters, m
"""
Nosecone_length = 1.775 #Nosecone length
Body_dia = 0.375 #body diameter
Body_len = 7.865 #body tube length
CoM = 6.10 #for Griffin

"""
Derived variables - Geometry m
"""
total_length = Body_len+Nosecone_length
fineness = total_length/Body_dia


"""
Material Properties
"""
G_alu = 24E9 #Shear Modulus Aluminium, Pa
G_alu_psi = 3.7e6 #Shear Modulus Aluminium, Psi
density_alu = 2710

"""
Simlation variables - Physical
"""
Roughness =3e-6 #body surface roughness - corresponding to a painted surface
N_fins = 4 #number of rocket fins
angle_attack = 0 #angle of attack for simulation
speed_sound = 295 #speed of sound