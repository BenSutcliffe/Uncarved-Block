from msilib.schema import Class
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import Drag
import Classes
import Parachute
import Stability

dt = 0.0001  # Change this to alter the timestep
iterations = 4000000  # Change this to alter how many timesteps are calculated (Should be in the order 1x10^6 timesteps)
droguecount = 0  # Variables for deploying times
maincount = 0

ratio = 0.6609  # Change this to alter the mass ratio (Total dry mass/ Total wet mass)
mass = [28]  # Change to alter the initial wet mass
densities = [1.225, 1.112, 1.007, 0.9093, 0.8194, 0.7364, 0.6601, 0.5900, 0.5258, 0.4671, 0.4135, 0.1948, 0.08891, 0.04008, 0.01841, 0.003996, 0.001027, 0.0003097, 0.00008283, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846]  # Air density data
altitudes = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000]
time = []

burntime = 6.3
for i in range(0, iterations + 100):
  time.append(dt * i)

time = np.array(time)

drag = np.zeros(iterations + 100)  # the +100 is to make sure all lists are longer than needed
height = np.zeros(iterations + 100)
velocity = np.zeros(iterations + 100)
acceleration = np.zeros(iterations + 100)
bfin_drag = np.zeros(iterations + 100)
dynamic_pressure = np.zeros(iterations + 100)

height[0] = 0 #20000
velocity[0] = 0 #824

for i in range(0, iterations + 100):
  mass.append((1 - ratio)*mass[0])
mass = np.array(mass)


for i in range(0, iterations):

  if height[i] >= 0 and height[i]<60000 and time[i] < 30:

    for j in range(0, len(altitudes)):
       if height[i] >= altitudes[j] and height[i] < altitudes[j+1]:
        density = densities[j] + ((densities[j+1] - densities[j]) * ((height[i] - altitudes[j])/(altitudes[j+1] - altitudes[j])))
        break

    if velocity[i] >= 0:
      Machno = velocity[i]/295
      Reynold = Drag.Reynolds(Classes.total_length, velocity[i])
      R_critical = Drag.R_crit(Classes.Roughness, Classes.total_length)
      Coeff_friction = Drag.friction_Re(Reynold, R_critical, Classes.Roughness, Classes.total_length)
      
      Coeff_friction_total = Drag.friction_vel(Coeff_friction, Machno, Classes.Bodyone.Arearef(), Classes.fineness, (np.pi *Classes.Bodyone.areat()), (Classes.Base.X_t*10**-3), Classes.Base.area(), Classes.MAC_base)
      Coeff_body_pressure = Drag.Body_pressure_drag(Machno, Classes.Bodyone.Arearef(), (Classes.Cone.L*10^(-2)))
      Coeff_base_drag = Drag.fin_pressure_drag(Machno, Classes.Leadingangle)

      Coefficient_of_drag = Coeff_friction_total + Coeff_base_drag[1] + Coeff_body_pressure + (Coeff_base_drag[0] * (4 * (Classes.Base.X_t*10**(-3)) * (Classes.Base.s*10**(-2)))/Classes.Bodyone.Arearef())
      drag[i] = Coefficient_of_drag * ((velocity[i]) ** 2) * 0.5 * density * Classes.Bodyone.Arearef()

    elif velocity[i] < 0:


      if height[i] > Parachute.drogueheight:
        drag[i] = -1 * Parachute.cddown * ((velocity[i]) ** 2) * 0.5 * density * Classes.Bodyone.Arearef()
      elif height[i] < Parachute.drogueheight and height[i] > Parachute.mainheight:
        if droguecount <= Parachute.droguedeploytime/dt:
          droguecount += 1
        drag[i] = -1 * Parachute.cddrogue * ((velocity[i]) ** 2) * 0.5 * density * (droguecount/(Parachute.droguedeploytime/dt)) * Parachute.adrogue

      elif height[i] < Parachute.mainheight:
        if maincount <= Parachute.maindeploytime/dt:
          maincount += 1
        drag[i] = -1 * Parachute.cdmain * ((velocity[i]) ** 2) * 0.5 * density * (maincount/(Parachute.maindeploytime/dt)) * Parachute.amain

    if float(dt*i) <= burntime:
      mass[i+1] = mass[i] - (float(1500)/94500000) * (ratio) * mass[0]
      acceleration[i+1] = (float(1500) - float(drag[i]))/float(mass[i]) - 9.81  # Gravity
    else:
        acceleration[i + 1] = (float(0) - float(drag[i])) / float(mass[i]) - 9.81

    velocity[i+1] = float(velocity[i]) + float(acceleration[i]) * dt
    height[i+1] = float(height[i]) + float(velocity[i]) * dt
    dynamic_pressure[i+1] = density * (float(velocity[i]**2))
    print(time[i], height[i], velocity[i], drag[i])

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