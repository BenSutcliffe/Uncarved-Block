from Drag import *
import Classes
from Classes import Fins, Body
from Thrust import *
import Global
from Stability import X_sf, X_N, C_N, CNalphaN_subs, CNalphaN_super

A_f = Classes.Base
A_len = Classes.total_length
Rough = Classes.Roughness
A_body = Classes.Bodyone
Transition_G = Classes.Transition
cone = Classes.Cone
P_dia = Classes.Panthera_dia

N_1 = Classes.N
alpha = 0 #Global.angle_attack

def Griffin_stability(Machno, Aquila_length, Aquila_f, n_cone, Aquila_Body, Panthera_f, Panthera_body, Trans_boat, alpha_f, N_f):
    beta = np.sqrt(np.abs(1-(Machno)**2))
    fin_pos_Aquila = (Aquila_length*10**(2)) - Aquila_f.C_r
    fin_pos_Panthera = (Aquila_length*10**(2)) + Panthera_body.h - Panthera_f.C_r

    if Machno <= 0.8:
        CNalpha_a = CNalphaN_subs(N_f, Aquila_f.s, Panthera_body.Arearef(), Aquila_f.area(), beta, Aquila_f.angle_skew())
        X_1_a = (Aquila_f.X_f() + fin_pos_Aquila) # Finds the centre of pressure for finset relative to the top of the rocket

        CNalpha_f = CNalphaN_subs(N_f, Panthera_f.s, Panthera_body.Arearef(), Panthera_f.area(), beta, Panthera_f.angle_skew())
        X_1_f = (Panthera_f.X_f() + fin_pos_Panthera) # Finds the centre of pressure for finset relative to the top of the rocket

    elif Machno < 1.4:
        CNalphasubs_a = CNalphaN_subs(N_f, Aquila_f.s, Panthera_body.Arearef(), Aquila_f.area(), 0.6, Aquila_f.angle_skew())
        CNalphasuper_a = CNalphaN_super(N_f, Panthera_body.Arearef(), Aquila_f.area(), 0.97979, alpha_f)
        CNalpha_a = CNalphasubs_a + ((Machno - 0.8)/0.6) * (CNalphasuper_a - CNalphasubs_a)
        X_1subs_a = (Aquila_f.X_f() + fin_pos_Aquila)
        X_1super_a = (X_sf(Aquila_f.Mac_f, Aquila_f.s, Aquila_f.area(), 0.97979) + fin_pos_Aquila)
        X_1_a = X_1subs_a + ((Machno - 0.8)/0.6) * (X_1super_a - X_1subs_a)

        CNalphasubs_f = CNalphaN_subs(N_f, Panthera_f.s, Panthera_body.Arearef(), Panthera_f.area(), 0.6, Panthera_f.angle_skew())
        CNalphasuper_f = CNalphaN_super(N_f, Panthera_body.Arearef(), Panthera_f.area(), 0.97979, alpha_f)
        CNalpha_f = CNalphasubs_f + ((Machno - 0.8)/0.6) * (CNalphasuper_f - CNalphasubs_f)
        X_1subs_f = (Panthera_f.X_f() + fin_pos_Panthera)
        X_1super_f = (X_sf(Panthera_f.Mac_f, Panthera_f.s, Panthera_f.area(), 0.97979) + fin_pos_Panthera)
        X_1_f = X_1subs_f + ((Machno - 0.8)/0.6) * (X_1super_f - X_1subs_f)

    else:
        CNalpha_a = CNalphaN_super(N_f, Panthera_body.Arearef(), Aquila_f.area(), beta, alpha_f)
        X_1_a = (X_sf(Aquila_f.Mac_f, Aquila_f.s, Aquila_f.area(), beta) + fin_pos_Aquila)

        CNalpha_f = CNalphaN_super(N_f, Panthera_body.Arearef(), Panthera_f.area(), beta, alpha_f)
        X_1_f = (X_sf(Panthera_f.Mac_f, Panthera_f.s, Panthera_f.area(), beta) + fin_pos_Panthera)


    X_3_a = X_N(Aquila_Body.n_1, Aquila_Body.h) # Finds the centre of pressure for body relative to the top of the rocket
    X_4_a = n_cone.X_f() # Finds the centre of pressure for the nosecone relative to the top of the rocket

    X_3_f = X_N((Aquila_length*10**(2)), Panthera_body.h) # Finds the centre of pressure for body relative to the top of the rocket

    C_1_a = CNalpha_a * Aquila_f.K() #Finds the normal force coefficent for a finset
    C_3_a = C_N(Aquila_Body.areat(), Aquila_Body.Arearef(), alpha_f) #Finds the normal force coefficent for the body
    C_4_a = 0.094 #n_cone.C[alpha_f] Not valid in combined case

    C_1_f = CNalpha_f * Panthera_f.K() #Finds the normal force coefficent for a finset
    C_3_f = C_N(Panthera_body.areat(), Panthera_body.Arearef(), alpha_f) #Finds the normal force coefficent for the body

    C_2_g = Trans_boat.CN_tr()
    X_2_g = Trans_boat.X_tr() + (Aquila_length*10**(2))

    #Finds the overall centre of pressure as a weighted mean of component COP and force coefficients, as from OpenRocket Documentation
    X = ((X_1_a * C_1_a) + (X_3_a * C_3_a) + (X_4_a * C_4_a) + (C_2_g*X_2_g) + (X_3_f * C_3_f)+ (X_1_f * C_1_f))/(C_1_a + C_3_a + C_4_a + C_3_f + C_1_f + C_2_g)
    Cal = (X - Classes.CoM2)/Panthera_body.d
    return Cal #Stores overall COP



P_C_r = 80
P_h = 35
P_C_t = 60
P_f_t = 10
MAC_base = 0
P_len = 750
Pan_f = Fins(P_C_r, P_h, P_C_t, P_f_t, MAC_base, 1)
Pan_body = Body(P_dia, P_len, 0)

print(Griffin_stability(0.3, A_len, A_f, cone, A_body, Pan_f, Pan_body, Transition_G, alpha, N_1))


