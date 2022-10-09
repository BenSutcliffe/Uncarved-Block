from site import venv
from Drag import *
import Classes
from Classes import Fins, Body
from stability import MAC, c_g, MAC_x, c_LE, flutter_NACA, G_alu_psi
from Global import *
from stability import X_sf, X_N, C_N, CNalphaN_subs, CNalphaN_super, pressure_position_transonic, C_N_force, F_fin_N
from scipy.optimize import differential_evolution

tot_len = Classes.total_length
Rough = Classes.Roughness
A_body = Classes.Bodyone
cone = Classes.Cone
P_dia = Classes.Body_dia
COM_est = Classes.CoM
body_len = Classes.Body_len

N_1 = Classes.N
alpha = 0 #Global.angle_attack

def Griffin_stability(Machno, total_length, n_cone, Panthera_f, Panthera_body, alpha_f, N_f, COM):
    beta = np.sqrt(np.abs(1-(Machno)**2))
    fin_pos_Panthera = (total_length*10**(2)) - Panthera_f.C_r

    if Machno <= 0.8:
        CNalpha_f = CNalphaN_subs(N_f, Panthera_f.s, Panthera_body.Arearef(), Panthera_f.area(), beta, Panthera_f.angle_skew())
        X_1_f = (Panthera_f.X_f() + fin_pos_Panthera) # Finds the centre of pressure for finset relative to the top of the rocket

    elif Machno < 1:
        CNalpha_f = CNalphaN_subs(N_f, Panthera_f.s, Panthera_body.Arearef(), Panthera_f.area(), beta, Panthera_f.angle_skew())
        X_1_f = pressure_position_transonic(MAC(Panthera_f, c_g), Panthera_f.s, Panthera_f.area(), Machno) + fin_pos_Panthera + MAC_x(Panthera_f, c_LE)
    
    elif Machno < 1.2:
        CNalpha_f = CNalphaN_super(N_f, Panthera_body.Arearef(), Panthera_f.area(), beta, alpha_f)
        X_1_f = pressure_position_transonic(MAC(Panthera_f, c_g), Panthera_f.s, Panthera_f.area(), Machno) + fin_pos_Panthera + MAC_x(Panthera_f, c_LE)

    else:
        CNalpha_f = CNalphaN_super(N_f, Panthera_body.Arearef(), Panthera_f.area(), beta, alpha_f)
        X_1_f = (X_sf(MAC(Panthera_f, c_g), Panthera_f.s, Panthera_f.area(), beta) + fin_pos_Panthera + MAC_x(Panthera_f, c_LE))

    X_4_a = n_cone.X_f() # Finds the centre of pressure for the nosecone relative to the top of the rocket

    X_3_f = X_N(n_cone.L, (Panthera_body.h)) # Finds the centre of pressure for body relative to the top of the rocket

    C_4_a = n_cone.CombinedC #n_cone.C[alpha_f] Not valid in combined case

    C_1_f = CNalpha_f * Panthera_f.K() #Finds the normal force coefficent for a finset
    C_3_f = C_N(Panthera_body.areat(), Panthera_body.Arearef(), alpha_f) #Finds the normal force coefficent for the body

    #Finds the overall centre of pressure as a weighted mean of component COP and force coefficients, as from OpenRocket Documentation
    X = ((X_4_a * C_4_a) + (X_3_f * C_3_f)+ (X_1_f * C_1_f))/(C_4_a + C_3_f + C_1_f)
    
    Cal = (X - COM)/Panthera_body.d
    return Cal #Stores overall COP

'''
density_alu = 2710
length = 786.5
h_s = 19
chord_r = 85
chord_t = 69
M_crit = 5.5
Aratio = 2*(h_s**2)/(0.5*(chord_r + chord_t)*h_s)
lambda_0 = chord_t/chord_r

Pant_f = Fins(chord_r, h_s, chord_t, (chord_r - chord_t), P_dia, 1)
c = MAC(Pant_f, c_g)/100
P_thick = flutter_NACA(1588.34, 295, 3235.974, 94321.68, c, Aratio, G_alu, lambda_0)
Pant_body = Body(P_dia, length, 0)
mass_t = P_thick * 0.5 * ((chord_r + chord_t)*0.01) * (h_s * 0.01) * density_alu *6 
COM_p = 610
print(P_thick)
'''

def objective(v):
    chord_r, chord_t, h_s = v
    length = body_len
    density_alu = 2710
    lambda_0 = chord_t/chord_r

    Aratio = (h_s**2)/(0.5*(chord_r + chord_t)*h_s)
    Pant_f = Fins(chord_r, h_s, chord_t, (chord_r - chord_t), P_dia, 1)

    c = MAC(Pant_f, c_g)/100
    P_thick = flutter_NACA(885, 295, 3235.974, 94321.68, c, Aratio, G_alu_psi, lambda_0)*1.5

    Pant_body = Body(P_dia, length, 0)

    mass_t = P_thick * 0.5 * ((chord_r + chord_t)*0.01) * (h_s * 0.01) * density_alu * N_1
    
    speeds = np.linspace(0.3, 5.5, 30)
    for velocity in speeds:
        M_stab = Griffin_stability(velocity, tot_len, cone, Pant_f, Pant_body, alpha, (N_1*0.913), COM_est)
        if M_stab < 1 or chord_t > (chord_r-15):
            output = 10e9
            return output
        else:
            pass
    output = mass_t
    return output
    

bounds = [[1, 150], [1, 150], [1, 70]]
result = differential_evolution(objective, bounds)
print('Status : %s' % result['message'])
print('Total Evaluations: %d' % result['nfev'])
# evaluate solution
solution = result['x']
evaluation = objective(solution)
print('Solution: f(%s) = %.5f' % (solution, evaluation))

chord_r = solution[0]
chord_t = solution[1]
h_s = solution[2]
alpha = 5
d_abso = alpha*np.pi/180
length = body_len
Pant_f = Fins(chord_r, h_s, chord_t, (chord_r - chord_t), P_dia, 1)
lambda_0 = chord_t/chord_r
c = MAC(Pant_f, c_g)/100
Aratio = (h_s**2)/(0.5*(chord_r + chord_t)*h_s)
Pant_body = Body(P_dia, length, 0)
P_thick = flutter_NACA(885, 295, 3235.974, 94321.68, c, Aratio, G_alu_psi, lambda_0)*1.5
CNalpha_f = CNalphaN_super(N_1, Pant_body.Arearef(), Pant_f.area(), 5.41, alpha)
C_N = C_N_force(CNalpha_f, d_abso)
F_fin = F_fin_N(C_N, max_q_rho, Pant_body.Arearef(), max_q_velo, d_abso)
print(P_thick, F_fin)