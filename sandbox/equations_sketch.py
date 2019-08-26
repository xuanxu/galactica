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
halo_radio_kpc = 110  # radio of the halo in kiloparsecs
disk_height_kpc = 0.2  # height of the disk in kiloparsecs
G = 0.44985  # Gravitational constant in Kpc^3/(10^9Msun * 10^7yrs)

def model():
    # Initial values of the system (y)
    initial_values = collections.defaultdict(float)

    gas_H = initial_values['gas_H']
    n = 1.5
    gas_H_n =  gas_H ** n
    cloud_H = initial_values['cloud_H']
    s1h = initial_values['s1h']
    s2h = initial_values['s2h']


    Kh1, Kh2 = star_formation_factor_halo()
    f = 1
    Wh = 0
    D1h = 0
    D2h = 0
    # Derivatives (ẏ)
    equations = {}

    equations['g_halo'] = -((Kh1 + Kh2) * gas_H_n) - (f * gas_H) + Wh
    equations['s_low_halo'] = 0.0
    equations['s_low_halo'] = (Kh1 * gas_H_n) - D1h
    equations['s_massive_halo'] = (Kh2 * gas_H_n) - D2h

    return equations


def star_formation_factor_halo():
    efficiency = 0.03  # epsilon_h from Ferrini et at, 1994, ApJ 427, 745
    factor = efficiency * (G / volume_halo())**0.5
    return [factor*0.5, factor*0.5]


def volume_halo(region_shape='square'):
    if region_shape == 'square':
        h = math.sqrt((halo_radio_kpc ** 2) - (region_galactocentric_radio_kpc ** 2))
        square_area = region_width_kpc * region_width_kpc
        return square_area * 2 * h
    elif region_shape == 'ring':
        h = math.sqrt((halo_radio_kpc ** 2) - (region_galactocentric_radio_kpc ** 2))
        half_ring_width = 0.5 * region_width_kpc
        ring_area = math.pi * (
            (region_galactocentric_radio_kpc + half_ring_width) ** 2 -
            (region_galactocentric_radio_kpc - half_ring_width) ** 2)
        return ring_area * 2 * h
    else:
        raise Exception("Wrong region shape. Allowed options: [square, ring]")
