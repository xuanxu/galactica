import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate import odeint


def initial_values():
  gas_H_disk = 1e11
  gas_H_halo = 1e12
  molecular_gas_disk = 0
  molecular_gas_halo = 0
  s_low_disk = 0
  s_massive_disk = 0
  s_low_halo = 0
  s_massive_halo = 0
  s_death_low_mass_disk = 0
  s_death_massive_disk = 0
  s_remnants_disk = 0
  Wd = 0
  y0 = [gas_H_disk, molecular_gas_disk, gas_H_halo, molecular_gas_halo]
  y0 = [gas_H_disk, molecular_gas_disk, s_death_low_mass_disk, s_death_massive_disk, s_remnants_disk]

  return y0


def disk_h(y, t, k):

  gas_D = y[0]# initial_values['gas_disk']
  c = y[1] # initial_values['molecular_gas_disk']
  gas_H = y[2] # initial_values['gas_halo']
  c_H = y[3] # initial_values['molecular_gas_halo']
  c2 = c ** 2

  n = 1.5
  gas_D_n = gas_D ** n

  Kc = 1 # star_formation_factor_cloud

  # star_formation_cloud_massive_stars_factor
  Ka1 = 1
  Ka2 = 1
  Ka_rest = 0.1

  # star_formation_cloud_collisions_factor
  Ks1 = 1
  Ks2 = 1
  Ks_rest = 0.1

  S2d = 1.0 # s_massive_disk dummy value
  Wd = 0

  # Death rates for low-mass (<8 Msun) and massive (>8 Msun) stars
  D1d = 0
  D2d = 0

  f = 1 # infall rate


  disk_H_gas = (-Kc * gas_D_n) + (Ka_rest * c * S2d) + (Ks_rest * c2) + (f * gas_H) + Wd
  disk_H_cloud = (Kc * gas_D_n) - ((Ka1 + Ka2 + Ka_rest) * c * S2d) - ((Ks1 + Ks2 + Ks_rest) * c2)
  disk_stars_low_mass = (Ks1 * c2) + (Ka1 * c * S2d) - D1d
  disk_stars_massive = (Ks2 * c2) + (Ka2 * c * S2d) - D2d
  disk_star_remnants = D1d + D2d - Wd

  dtdy = [disk_H_gas, disk_H_cloud, gas_H, c_H]
  dtdy = [disk_H_gas, disk_H_cloud, disk_stars_low_mass, disk_stars_massive, disk_star_remnants]
  return dtdy

# Define time points
t = np.linspace(0, 14, 10000)

# Define initial values
y0 = initial_values()

# Define parameters
params = 1

# Solve the system of differential equations
solution = odeint(disk_h, y0, t, args=(params,))

# Plot the results
plt.figure()
plt.plot(t, solution[:, 0], label='Gas Disk')
plt.plot(t, solution[:, 1], label='Molecular Gas Disk')
#plt.plot(t, solution[:, 2], label='Gas Halo')
#plt.plot(t, solution[:, 3], label='Molecular Gas Halo')
plt.plot(t, solution[:, 2], label='Low mass stars')
plt.plot(t, solution[:, 3], label='Massive stars')
plt.plot(t, solution[:, 4], label='Remnants')

plt.xlabel('Time')
plt.ylabel('Value in Msun')
plt.legend()
plt.show()


