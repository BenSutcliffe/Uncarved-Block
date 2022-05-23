from msilib.schema import Class
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import Drag
import Classes
import Thrust
import Global
from Stability import X_sf, X_N, C_N, CNalphaN_subs, CNalphaN_trans, CNalphaN_super

dt = 0.0001  # Change this to alter the timestep
iterations = 4000000  # Change this to alter how many timesteps are calculated (Should be in the order 1x10^6 timesteps)
droguecount = 0  # Variables for deploying times
maincount = 0

ratio = 0.6609  # Change this to alter the mass ratio (Total dry mass/ Total wet mass)
densities = [1.225, 1.112, 1.007, 0.9093, 0.8194, 0.7364, 0.6601, 0.5900, 0.5258, 0.4671, 0.4135, 0.1948, 0.08891, 0.04008, 0.01841, 0.003996, 0.001027, 0.0003097, 0.00008283, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846]  # Air density data
altitudes = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000]
time = []

burntime = 3.6
for i in range(0, iterations + 100):
  time.append(dt * i)

time = np.array(time)

drag = np.zeros(iterations + 100)  # the +100 is to make sure all lists are longer than needed
height = np.zeros(iterations + 100)
velocity = np.zeros(iterations + 100)
acceleration = np.zeros(iterations + 100)
bfin_drag = np.zeros(iterations + 100)
dynamic_pressure = np.zeros(iterations + 100)
mass = np.zeros(iterations + 100)
calibre = np.zeros(iterations + 100)

height[0] = 0 #20000
velocity[0] = 0 #824
mass[0] = Global.start_mass



for i in range(0, iterations):

  if height[i] >= 0 and height[i]<150000 and time[i] < 80:

    for j in range(0, len(altitudes)):
       if height[i] >= altitudes[j] and height[i] < altitudes[j+1]:
        density = densities[j] + ((densities[j+1] - densities[j]) * ((height[i] - altitudes[j])/(altitudes[j+1] - altitudes[j])))
        break

    if velocity[i] >= 0:
      if height[i] <= 35000:
        a_speedsound = 340 - ((height[i]/35000)*45)
      else:
        a_speedsound = 295

      Machno = velocity[i]/a_speedsound


      beta = np.sqrt(np.abs(1-(Machno)**2))
      Reynold = Drag.Reynolds(Classes.total_length, velocity[i])
      R_critical = Drag.R_crit(Classes.Roughness, Classes.total_length)
      Coeff_friction = Drag.friction_Re(Reynold, R_critical, Classes.Roughness, Classes.total_length)
      
      Coeff_friction_total = Drag.friction_vel(Coeff_friction, Machno, Classes.Bodyone.Arearef(), Classes.fineness, (np.pi *Classes.Bodyone.areat()), (Classes.Base.X_t*10**-3), Classes.Base.area(), Classes.MAC_base)
      Coeff_body_pressure = Drag.Body_pressure_drag(Machno, Classes.Bodyone.Arearef(), (Classes.Cone.L*10**(-2)))
      Coeff_base_drag = Drag.fin_pressure_drag(Machno, Classes.Leadingangle)

      Coefficient_of_drag = Coeff_friction_total + Coeff_base_drag[1] + Coeff_body_pressure + (Coeff_base_drag[0] * (4 * (Classes.Base.X_t*10**(-3)) * (Classes.Base.s*10**(-2)))/Classes.Bodyone.Arearef())
      drag[i] = Coefficient_of_drag * ((velocity[i]) ** 2) * 0.5 * density * Classes.Bodyone.Arearef()

      if Machno <= 0.8:
        CNalpha = CNalphaN_subs(Classes.N, Classes.Base.s, Classes.Bodyone.Arearef(), Classes.Base.area(), beta, Classes.Base.gamma)
        X_1 = (Classes.Base.X_f() + 224) # Finds the centre of pressure for finset relative to the top of the rocket
      elif Machno <= 1.2:
        CNalpha = CNalphaN_trans(Classes.N, Classes.Bodyone.Arearef(), Classes.Base.area(), beta, Global.angle_attack)
        X_1 = (Classes.Base.X_f() + 224) # Finds the centre of pressure for finset relative to the top of the rocket
      else:
        CNalpha = CNalphaN_super(Classes.N, Classes.Bodyone.Arearef(), Classes.Base.area(), beta, Global.angle_attack)
        X_1 = (X_sf(27.282, Classes.Base.s, Classes.Base.area(), beta) + 224)

      X_3 = X_N(Classes.Bodyone.n_1, Classes.Bodyone.h) # Finds the centre of pressure for body relative to the top of the rocket
      X_4 = Classes.Cone.X_f() # Finds the centre of pressure for the nosecone relative to the top of the rocket

      C_1 = CNalpha * Classes.Base.K() #Finds the normal force coefficent for a finset
      C_3 = C_N(Classes.Bodyone.areat(), Classes.Bodyone.Arearef(), Global.angle_attack) #Finds the normal force coefficent for the body
      C_4 = Classes.Cone.C[Global.angle_attack] #Finds the normal force coefficent for the nosecone

      #Finds the overall centre of pressure as a weighted mean of component COP and force coefficients, as from OpenRocket Documentation
      X = ((X_1 * C_1) + (X_3 * C_3) + (X_4 * C_4))/(C_1 + C_3 + C_4)

      Cal = (X - Classes.CoM)/Classes.Bodyone.d
      calibre[i] = Cal #Stores overall COP

    elif velocity[i] < 0:
      if height[i] > Global.drogueheight:
        drag[i] = -1 * Global.cddown * ((velocity[i]) ** 2) * 0.5 * density * Classes.Bodyone.Arearef()
      elif height[i] < Global.drogueheight and height[i] > Global.mainheight:
        if droguecount <= Global.droguedeploytime/dt:
          droguecount += 1
        drag[i] = -1 * Global.cddrogue * ((velocity[i]) ** 2) * 0.5 * density * (droguecount/(Global.droguedeploytime/dt)) * Global.adrogue

      elif height[i] < Global.mainheight:
        if maincount <= Global.maindeploytime/dt:
          maincount += 1
        drag[i] = -1 * Global.cdmain * ((velocity[i]) ** 2) * 0.5 * density * (maincount/(Global.maindeploytime/dt)) * Global.amain

    if float(dt*i) <= burntime:
      thrust = Thrust.Thrust_curve(time[i])
      mass[i+1] = mass[i] - ((Global.start_mass-Global.end_mass)*dt/burntime)
      acceleration[i+1] = (float(thrust) - float(drag[i]))/float(mass[i]) - 9.81  # Gravity
    else:
        acceleration[i + 1] = (float(0) - float(drag[i]))/Global.end_mass - 9.81

    velocity[i+1] = float(velocity[i]) + float(acceleration[i]) * dt
    height[i+1] = float(height[i]) + float(velocity[i]) * dt
    dynamic_pressure[i+1] = density * (float(velocity[i]**2))
    print(time[i], height[i], velocity[i], drag[i], mass[i])

maxheight = np.sort(height)
print("Maximum height: ", maxheight[-1])
maxvelocity = np.sort(velocity)
print("Maximum velocity: ", maxvelocity[-1])
maxaccel = np.sort(acceleration)
print("Maximum acceleration: ", maxaccel[-1])
maxdrag = np.sort(drag)
print("Maximum drag: ", maxdrag[-1])
maxq = np.sort(dynamic_pressure)
print("")
print("Maximum Dynamic Pressure:", maxq[-1])

maxq = 0
for i in range(len(dynamic_pressure)):
    if float(dynamic_pressure[i]) > float(dynamic_pressure[maxq]):
        maxq = i

print("/nMax Q occurs at an altitude of", height[maxq], "and a speed of", velocity[maxq])

plt.plot(time[0:len(height):10], height[0:len(height):10])
plt.legend("height")
plt.show()

plt.plot(time[0:len(height):10], velocity[0:len(height):10])
plt.legend("velocity")
plt.show()

plt.plot(time[0:len(height):10], calibre[0:len(height):10])
plt.legend("Stability")
plt.show()