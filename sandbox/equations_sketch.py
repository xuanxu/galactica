import math
import collections

'''
  gas_H = diffuse gas in the galactic Halo
  gas_D = diffuse gas in the galactic Disc
  cloud_H = molecular gas in the halo
  cloud_D = molecular gas in the disc
  f = infall rate

  Stars are divided in two groups:
    1) Low/intermediate mass stars (m <= 4 Msun)
    2) Massive stars (m > 4 Msun)
  So the 1 and 2 subscripts refer to this groups.

  The d and h subscripts correspond to disc and halo.

  Kh = proportionality factor of the SF in the halo
  Kc = proportionality factor of the SF in the cloud formation
  Ks = proportionality factor of the SF in the cloud-cloud collision
  Ka = proportionality factor of the SF in the cloud-massive stars interactions
  Since stars are divided in two groups, the parameters
  involving Star Formation are divided in two groups too:
  Kh = Kh1 + Kh2
  Kc = Kc1 + Kc2
  Ks = Ks1 + Ks2 + Ks'
  Ka = Ka1 + Ka2 + Ka'
     where Ks' and Ka' refer to the restitution of diffuse gas
     due to the collision and interaction processes

  D1 = death rates for low/intermediate mass stars
  D2 = death rates for massive stars

  Wd = Restitution rate in the disc
  Wh = Restitution rate in the halo
'''

# Constants
region_width_kpc = 1  # width of the region in kiloparsecs
region_galactocentric_radio_kpc = 8  # distance of the region to the galactic center in kiloparsecs
halo_radio_kpc = 260.996  # radio of the halo in kiloparsecs
disk_height_kpc = 0.2  # height of the disk in kiloparsecs
G = 0.44985  # Gravitational constant in Kpc^3/(10^9Msun * 10^7yrs)
virial_mass = 1e12  # Virial mass of the dark matter halo
disc_barionic_mass = 1e11  # Mass of the barionic disc


def model():
    # Initial values of the system (y)
    initial_values = collections.defaultdict(float)
    initial_values['gas_halo'] = virial_mass
    initial_values['gas_disk'] = disc_barionic_mass

    gas_H = initial_values['gas_halo']
    gas_D = initial_values['gas_disk']
    n = 1.5
    gas_H_n = gas_H ** n
    gas_D_n = gas_D ** n
    molecular_gas_H = initial_values['molecular_gas_halo']
    molecular_gas_D = initial_values['molecular_gas_disk']
    S1h = initial_values['s_low_halo']
    S2h = initial_values['s_massive_halo']
    S1d = initial_values['s_low_disk']
    S2d = initial_values['s_massive_disk']

    Kh1, Kh2 = star_formation_factor_halo()
    Kc = star_formation_factor_cloud()
    Ka1, Ka2, Ka_rest = star_formation_cloud_massive_stars_factor()
    Ks1, Ks2, Ks_rest = star_formation_cloud_collisions_factor()
    f = 1
    Wd = 0
    Wh = 0
    D1d = 0
    D2d = 0
    D1h = 0
    D2h = 0
    c = molecular_gas_D
    c2= c ** 2

    # Derivatives (ẏ)
    equations = {'halo': {}, 'disk': {}}


    equations['disk']['gas'] = (-Kc * gas_D_n) + (Ka_rest * c * S2d) + (Ks_rest * c2) + (f * gas_H) + Wd
    equations['disk']['cloud'] = (Kc * gas_D_n) - ((Ka1 + Ka2 + Ka_rest) * c * S2d) - ((Ks1 + Ks2 + Ks_rest) * c2)
    equations['disk']['stars_low'] = (Ks1 * c2) + (Ka1 * c * S2d) - D1d
    equations['disk']['stars_massive'] = (Ks2 * c2) + (Ka2 * c * S2d) - D2d
    equations['disk']['remnants'] = D1d + D2d - Wd
    equations['halo']['gas'] = -((Kh1 + Kh2) * gas_H_n) - (f * gas_H) + Wh
    equations['halo']['cloud'] = 0.0
    equations['halo']['stars_low'] = (Kh1 * gas_H_n) - D1h
    equations['halo']['stars_massive'] = (Kh2 * gas_H_n) - D2h
    equations['halo']['remnants'] = D1h + D2h - Wh

    return equations


def integrator(state):
  pass


def star_formation_factor_halo():
    efficiency = 2.173  # epsilon_h computed for a best value K_h = 9e-3 able to reproduce SFR and abundances of MWG halo
    factor = efficiency * (G / volume_halo())**0.5
    return [factor*0.5, factor*0.5]


def star_formation_factor_cloud():
  return 1


def star_formation_cloud_massive_stars_factor():
  return [1, 1, 0.1]


def star_formation_cloud_collisions_factor():
  return [1, 1, 0.1]


def volume_halo(region_shape='square'):
    h = math.sqrt((halo_radio_kpc ** 2) - (region_galactocentric_radio_kpc ** 2))
    if region_shape == 'square':
        square_area = region_width_kpc * region_width_kpc
        return square_area * 2 * h
    elif region_shape == 'ring':
        half_ring_width = 0.5 * region_width_kpc
        ring_area = math.pi * (
            (region_galactocentric_radio_kpc + half_ring_width) ** 2 -
            (region_galactocentric_radio_kpc - half_ring_width) ** 2)
        return ring_area * 2 * h
    else:
        raise Exception("Wrong region shape. Allowed options: [square, ring]")


def volume_disk(region_shape='square'):
    if region_shape == 'square':
        square_area = region_width_kpc * region_width_kpc
        return square_area * disk_height_kpc
    elif region_shape == 'ring':
        half_ring_width = 0.5 * region_width_kpc
        ring_area = math.pi * (
            (region_galactocentric_radio_kpc + half_ring_width) ** 2 -
            (region_galactocentric_radio_kpc - half_ring_width) ** 2)
        return ring_area * disk_height_kpc
    else:
        raise Exception("Wrong region shape. Allowed options: [square, ring]")
