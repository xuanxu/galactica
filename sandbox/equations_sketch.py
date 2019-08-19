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
radio_kpc = 8  # radio of the region in kiloparsecs
G = 0.44985  # Gravitational constant in Kpc^3/(10^9Msun * 10^7yrs)


# Initial values of the system (y)
initial_values = []

gas_H = initial_values[0]
cloud_H = initial_values[1]
s1h = initial_values[2]
s2h = initial_values[3]


# Derivatives (ẏ)
equations = []

equations[0] = -((Kh1 + Kh2) * gas_H_n) - (f * gas_H) + Wh
equations[1] = 0.0
equations[2] = (Kh1 * gas_H_n) - D1h
equations[3] = (Kh2 * gas_H_n) - D2h


def gas_H_n:
    n = 1.5
    return gas_H ** n


def Kh():
    eps_h * (G / Vh)**0.5
