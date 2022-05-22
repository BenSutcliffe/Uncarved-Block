import numpy as np

from Classes import Nosecone_length

print(Nosecone_length)

import matplotlib.pyplot as plt
import numpy as np

dt = 0.0001  # Change this to alter the timestep
iterations = 4000000  # Change this to alter how many timesteps are calculated (Should be in the order 1x10^6 timesteps)

cddown = 0.4  # Change this to alter the Cd on descent
adown = 0.01539  # Change this to alter the effective area on descent
Length = np.arange(241)
A_Fin_base = 0.052
A_Fin_top = 2.4e-3
fin_thickness_base = 4e-3
MAC_base = 27.282e-2
Leadingangle_base = 26.6
Leadingangle_top = 33.7
Roughness =3e-6
totatlen = 2.4
noseconelen = 40e-2
maximumdiameter = 13e-2
A_ref = np.pi * (maximumdiameter/2)**2
A_body_wet = 2*np.pi*(maximumdiameter/2)*totatlen
A_fin_wet_base = 8 * A_Fin_base
A_fin_wet_top = 8 * A_Fin_top


cddrogue = 0.5  # Change this to alter the Cd of the drogue
adrogue = 0.0818 # Change this to alter the area of the drogue
amain = 1 # Change this to alter the area of the main chute
cdmain = 1.2  # Change this to alter the area of the main chute
drogueheight = 10000  # Change this to alter the drogue deploy height
mainheight = 1000  # Change this to alter the main chute deploy height
droguedeploytime = 10  # Change this to alter the drogue deploy time
maindeploytime = 10  # Change this to alter the main chute deploy time
ratio = 0.6609  # Change this to alter the mass ratio (Total dry mass/ Total wet mass)
mass = [24]  # Change to alter the initial wet mass
densities = [1.225, 1.112, 1.007, 0.9093, 0.8194, 0.7364, 0.6601, 0.5900, 0.5258, 0.4671, 0.4135, 0.1948, 0.08891, 0.04008, 0.01841, 0.003996, 0.001027, 0.0003097, 0.00008283, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846, 0.00001846]  # Air density data
altitudes = [0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000, 25000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000, 110000, 120000, 130000, 140000, 150000]

time = []
droguecount = 0  # Variables for deploying times
maincount = 0


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

def Reynolds(length, velocity, viscosity = 1.5e-5):
    R = (length * velocity)/viscosity
    return R

def R_crit(Roughness, Length):
    R_crit = 51 * (Roughness/Length)**(-1.039)
    return R_crit


for i in range(0, iterations + 100):
  mass.append((1 - ratio)*mass[0])
mass = np.array(mass)
#force = np.array(force)

for i in range(0, iterations):

  if height[i] >= 0 and height[i]<60000 and time[i] < 30:

    for j in range(0, len(altitudes)):
       if height[i] >= altitudes[j] and height[i] < altitudes[j+1]:
        density = densities[j] + ((densities[j+1] - densities[j]) * ((height[i] - altitudes[j])/(altitudes[j+1] - altitudes[j])))
        break

    if velocity[i] >= 0:
      #drag[i] = cdup * ((velocity[i]) ** 2) * 0.5 * density * aup
      Reynold = Reynolds(totatlen, velocity[i])
      R_critical = R_crit(Roughness, totatlen)
      Coeff_friction = 0
      if Reynold <= 1e4:
          Coeff_friction = 1.48e-2
      elif Reynold > 1e4 and Reynold < R_critical:
          Coeff_friction = ((1.5*np.log(Reynold))-5.6)**(-2)
      elif Reynold >= R_critical:
          Coeff_friction = 0.032 * ((Roughness/totatlen)**0.2)


      Machno = velocity[i]/295
      Skin_friction_coeff = 0
      if Machno < 0.8:
          Skin_friction_coeff = Coeff_friction * (1 - ((0.1) * (Machno ** 2)))
      elif Machno > 0.8:
          Supersonic_skin_friction_1 = Coeff_friction * (1 + ((0.15) * (Machno ** 2)))**(-0.58)
          Supersonic_skin_friction_2 = Coeff_friction * (1 + ((0.18) * (Machno ** 2)))**(-1)
          if Supersonic_skin_friction_1 > Supersonic_skin_friction_2:
              Skin_friction_coeff = Supersonic_skin_friction_1
          else:
              Skin_friction_coeff = Supersonic_skin_friction_2
      Coeff_friction_total = (Skin_friction_coeff/A_ref) * (((1 + (1/(2 * fineness))) * (A_body_wet)) + ((1 + (2 * fin_thickness_base/MAC_base)) * (A_fin_wet_base))) #+ ((1 + (2 * fin_thickness_top/MAC_top)) * (A_fin_wet_top)))

      #D_friction = Coeff_friction_total * ((velocity[i]) ** 2) * 0.5 * density * (A_ref)

      Coeff_body_pressure = 0
      if Machno < 0.8:
          Coeff_body_pressure = 0
      elif Machno < 1.3 and Machno >= 0.8:
          eta = np.arctan(maximumdiameter/(2*noseconelen))
          Coeff_body_pressure = np.sin(eta)
      elif Machno >= 1.3:
          eta = np.arctan(maximumdiameter/(2*noseconelen))
          Coeff_body_pressure = (2.1 * ((np.sin(eta))**2)) + ((np.sin(eta))/(2 * np.sqrt((Machno**2) - 1)))

      #D_body_pressure = Coeff_body_pressure * ((velocity[i]) ** 2) * 0.5 * density * A_ref

      Coeff_fin_pressure_LE = 0  #Assuming a rounded leading edge
      if Machno < 0.9:
          Coeff_fin_pressure_LE = ((1 - (Machno**2))**(-0.417)) - 1
      elif Machno >= 0.9 and Machno <1:
          Coeff_fin_pressure_LE = 1 - (1.785 * (Machno-0.9))
      elif Machno >= 1:
          Coeff_fin_pressure_LE = 1.214 - (0.502/(Machno**2)) + (0.1095/(Machno**4))

      Coeff_fin_pressure_LE_base = Coeff_fin_pressure_LE * (np.cos(Leadingangle_base * np.pi/180))**2
      Coeff_fin_pressure_LE_top = Coeff_fin_pressure_LE * (np.cos(Leadingangle_top * np.pi / 180))**2

      Coeff_fin_pressure_TE = 0
      if Machno < 1:
          Coeff_fin_pressure_TE = 0.5*(0.12 + (0.13 * (Machno**2)))
      elif Machno >= 1:
          Coeff_fin_pressure_TE = 0.5*(0.25/Machno)

      Coeff_fin_pressure_base = Coeff_fin_pressure_LE_base
      Coeff_fin_pressure_top = 0 #Coeff_fin_pressure_LE_top

      D_fin_pressure_base = Coeff_fin_pressure_base * ((velocity[i]) ** 2) * 0.5 * density * (4 * fin_thickness_base * 20e-2)
      #D_fin_pressure_top = Coeff_fin_pressure_top * ((velocity[i]) ** 2) * 0.5 * density * (4 * fin_thickness_top * 3e-2)


      Coeff_base_drag = Coeff_fin_pressure_TE
      #D_base = Coeff_base_drag * ((velocity[i]) ** 2) * 0.5 * density * A_ref

      #print("Friction: {}, Body pressure: {}, Fin base pressure: {}, Fin Top Pressure: {}, Base Drag: {}, Velocity: {}".format(D_friction, D_body_pressure, D_fin_pressure_base, D_fin_pressure_top, D_base, velocity[i]))

      Coefficient_of_drag = Coeff_friction_total + (Coeff_base_drag * A_ref/A_ref) + (Coeff_body_pressure * A_ref/A_ref) + (Coeff_fin_pressure_base * (4 * fin_thickness_base * 20e-2)/A_ref) + (Coeff_fin_pressure_top * (4 * fin_thickness_top * 3e-2)/A_ref)

      #drag[i] = D_friction + D_body_pressure + D_fin_pressure_base + D_fin_pressure_top + D_base
      drag[i] = Coefficient_of_drag * ((velocity[i]) ** 2) * 0.5 * density * A_ref
      bfin_drag[i] = (D_fin_pressure_base/4)

    elif velocity[i] < 0:


      if height[i] > drogueheight:
        drag[i] = -1 * cddown * ((velocity[i]) ** 2) * 0.5 * density * adown
      elif height[i] < drogueheight and height[i] > mainheight:
        if droguecount <= droguedeploytime/dt:
          droguecount += 1
        drag[i] = -1 * cddrogue * ((velocity[i]) ** 2) * 0.5 * density * (droguecount/(droguedeploytime/dt)) * adrogue

      elif height[i] < mainheight:
        if maincount <= maindeploytime/dt:
          maincount += 1
        drag[i] = -1 * cdmain * ((velocity[i]) ** 2) * 0.5 * density * (maincount/(maindeploytime/dt)) * amain

    if float(dt*i) <= burntime:
      mass[i+1] = mass[i] - (float(1500)/94500000) * (ratio) * mass[0]
      acceleration[i+1] = (float(1500) - float(drag[i]))/float(mass[i]) - 9.81  # Gravity
    else:
        acceleration[i + 1] = (float(0) - float(drag[i])) / float(mass[i]) - 9.81

    velocity[i+1] = float(velocity[i]) + float(acceleration[i]) * dt
    height[i+1] = float(height[i]) + float(velocity[i]) * dt
    dynamic_pressure[i+1] = density * (float(velocity[i]**2))
    print(time[i], height[i], velocity[i], drag[i])

plt.plot(time[0:len(height):10], height[0:len(height):10])
plt.legend("height")
plt.show()

plt.plot(time[0:len(height):10], velocity[0:len(height):10])
plt.legend("velocity")
plt.show()

plt.plot(time[0:len(height):10], acceleration[0:len(height):10])
plt.legend("acceleration")
plt.show()

plt.plot(time[0:len(height):10], drag[0:len(height):10])
plt.legend("drag")
plt.show()

plt.plot(time[0:len(height):10], mass[0:len(height):10])
plt.legend("mass")
plt.show()

plt.plot(time[0:len(height):10], dynamic_pressure[0:len(height):10])
plt.legend("q dynamic pressure")
plt.show()

plt.plot(height[0:9000:10], velocity[0:9000:10])
plt.legend("Velocity vs Height")
plt.show()

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