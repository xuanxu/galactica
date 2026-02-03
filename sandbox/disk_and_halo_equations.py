import collections
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def settings():
  s = collections.defaultdict(float)
  s["region_width_kpc"] = 1  # width of the region in kiloparsecs
  s["region_galactocentric_radio_kpc"] = 8  # distance of the region to the galactic center in kiloparsecs
  s["halo_radio_kpc"] = 260.996  # radio of the halo in kiloparsecs
  s["disk_height_kpc"] = 0.2 # height of the galactic disk in kiloparsecs
  s["G"] = 0.44985  # Gravitational constant in Kpc^3/(10^9Msun * 10^7yrs)

  return s

def initial_values(radius=1):
  # Disk
  disk_H_gas = 1e11
  disk_H2_gas = 0.0
  disk_low_mass_stars = 0.0
  disk_massive_stars = 0.0
  disk_stars_remnants = 0.0
  disk_death_low_mass_stars = 0.0
  disk_death_massive_stars = 0.0
  # Halo
  halo_H_gas = 1e12
  halo_H2_gas = 0.0
  halo_low_mass_stars = 0.0
  halo_massive_stars = 0.0
  halo_stars_remnants = 0.0

  return [disk_H_gas, disk_H2_gas, disk_low_mass_stars, disk_massive_stars, disk_stars_remnants,
          halo_H_gas, halo_low_mass_stars, halo_massive_stars, halo_stars_remnants]


def cloud_formation_factor():
  return 1


def star_formation_cloud_massive_stars_collisions_factors():
  low_stars_factor = 1.0
  massive_stars_factor = 1.0
  restitution_to_diffuse_gas = 0.1
  return [low_stars_factor, massive_stars_factor, restitution_to_diffuse_gas]


def star_formation_cloud_cloud_collisions_factors():
  low_stars_factor = 1.0
  massive_stars_factor = 1.0
  restitution_to_diffuse_gas = 0.1
  return [low_stars_factor, massive_stars_factor, restitution_to_diffuse_gas]


def star_formation_in_halo_factors():
  s = settings()
  efficiency = 2.173  # epsilon_h computed for a best value K_h = 9e-3 able to reproduce SFR and abundances of MWG halo
  factor = efficiency * (s["G"] / volume_halo())**0.5

  low_stars_factor = factor*0.5
  massive_stars_factor = factor*0.5
  return [low_stars_factor, massive_stars_factor]


def volume_halo(region_shape='ring'):
  s = settings()
  h = np.sqrt((s["halo_radio_kpc"] ** 2) - (s["region_galactocentric_radio_kpc"] ** 2))
  if region_shape == 'square':
      square_area = s["region_width_kpc"] ** 2
      return square_area * 2 * h
  elif region_shape == 'ring':
      half_ring_width = 0.5 * s["region_width_kpc"]
      ring_area = np.pi * (
          (s["region_galactocentric_radio_kpc"] + half_ring_width) ** 2 -
          (s["region_galactocentric_radio_kpc"] - half_ring_width) ** 2)
      return ring_area * 2 * h
  else:
      raise Exception("Wrong region shape. Allowed options: [square, ring]")


def volume_disk(region_shape='ring'):
  s = settings()
  if region_shape == 'square':
      square_area = s["region_width_kpc"] ** 2
      return square_area * s["disk_height_kpc"]
  elif region_shape == 'ring':
      half_ring_width = 0.5 * s["region_width_kpc"]
      ring_area = np.pi * (
          (s["region_galactocentric_radio_kpc"] + half_ring_width) ** 2 -
          (s["region_galactocentric_radio_kpc"] - half_ring_width) ** 2)
      return ring_area * s["disk_height_kpc"]
  else:
      raise Exception("Wrong region shape. Allowed options: [square, ring]")


def parameters(radius=1):
  params = collections.defaultdict(float)
  params["radius"] = radius
  params["infall_rate"] = 0.1
  params["cloud_formation_factor"] = cloud_formation_factor()
  params["star_formation_in_halo"] = star_formation_in_halo_factors() # Kh1 + Kh2
  params["star_formation_cloud_cloud_collisions"] = star_formation_cloud_cloud_collisions_factors() # Ks1 + Ks2 + Ks_restitution_of_diffuse_gas
  params["star_formation_cloud_massive_stars_collisions"] = star_formation_cloud_massive_stars_collisions_factors() # Ka1 + Ka2 + Ka_restitution_of_diffuse_gas

  return params


def galaxy_system(t, values, params):

  disk_H_gas, disk_H2_gas, disk_low_mass_stars, disk_massive_stars, disk_stars_remnants, \
  halo_H_gas, halo_low_mass_stars, halo_massive_stars, halo_stars_remnants = values

  radius = params["radius"]

  n = 1.5
  disk_H_gas_n = disk_H_gas ** n
  disk_H2_gas_2 = disk_H2_gas ** 2
  halo_H_gas_n = halo_H_gas ** n

  Kc = params["cloud_formation_factor"]

  # SF factor by cloud - massive stars interaction
  Ka1 = params["star_formation_cloud_massive_stars_collisions"][0]
  Ka2 = params["star_formation_cloud_massive_stars_collisions"][1]
  Ka_rest = params["star_formation_cloud_massive_stars_collisions"][2]

  # SF factor by cloud - cloud collisions
  Ks1 = params["star_formation_cloud_cloud_collisions"][0]
  Ks2 = params["star_formation_cloud_cloud_collisions"][1]
  Ks_rest = params["star_formation_cloud_cloud_collisions"][2]

  Kh1 = params["star_formation_in_halo"][0]
  Kh2 = params["star_formation_in_halo"][1]

  S2d = 1.0 # s_massive_disk dummy value
  Wd = 0
  Wh = 0

  # Death rates for low-mass (<8 Msun) and massive (>8 Msun) stars
  D1d = 0
  D2d = 0
  D1h = 0
  D2h = 0


  infall_rate = params["infall_rate"]

  dt_disk_H_gas = (-Kc * disk_H_gas_n) + (Ka_rest * disk_H2_gas * S2d) + (Ks_rest * disk_H2_gas_2) + (infall_rate * halo_H_gas) + Wd
  dt_disk_H2_gas = (Kc * disk_H_gas_n) - ((Ka1 + Ka2 + Ka_rest) * disk_H2_gas * S2d) - ((Ks1 + Ks2 + Ks_rest) * disk_H2_gas_2)
  dt_disk_low_mass_stars = (Ks1 * disk_H2_gas_2) + (Ka1 * disk_H2_gas * S2d) - D1d
  dt_disk_massive_stars = (Ks2 * disk_H2_gas_2) + (Ka2 * disk_H2_gas * S2d) - D2d
  dt_disk_stars_remnants = D1d + D2d - Wd

  dt_halo_H_gas = -((Kh1 + Kh2) * halo_H_gas_n) - (infall_rate * halo_H_gas) + Wh
  dt_halo_low_mass_stars = (Kh1 * halo_H_gas_n) - D1h
  dt_halo_massive_stars = (Kh2 * halo_H_gas_n) - D2h
  dt_halo_stars_remnants = D1h + D2h - Wh

  dtdy = [ dt_disk_H_gas,
           dt_disk_H2_gas,
           dt_disk_low_mass_stars,
           dt_disk_massive_stars,
           dt_disk_stars_remnants,
           dt_halo_H_gas,
           dt_halo_low_mass_stars,
           dt_halo_massive_stars,
           dt_halo_stars_remnants ]

  return dtdy




time_span = [0, 13.2]
time_eval = np.linspace(time_span[0], time_span[1], 100)
radius = 1

solution = solve_ivp(galaxy_system, time_span, initial_values(), args=(parameters(radius),), t_eval=time_eval)

t = solution.t
y = solution.y[0]




plt.title("Galaxy evolution")
plt.plot(t, solution.y[0], label='Diffuse gas in disk')
plt.plot(t, solution.y[1], label='Molecular gas in disk')
plt.plot(t, solution.y[2], label='Low mass stars in disk')
plt.plot(t, solution.y[3], label='Massive stars in disk')
plt.plot(t, solution.y[4], label='Remnants in disk')
plt.plot(t, solution.y[5], label='Diffuse gas in halo')
plt.plot(t, solution.y[6], label='Low mass stars in halo')
plt.plot(t, solution.y[7], label='Massive stars in halo')
plt.plot(t, solution.y[8], label='Remnants in halo')
plt.xlabel("Time")
plt.ylabel("Msun")
plt.legend()
plt.grid()
plt.show()


