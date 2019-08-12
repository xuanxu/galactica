'''
  H = hydrogen
  f = infall rate
  gas_H = diffuse gas in the galactic Halo
  gas_D = diffuse gas in the galactic Disc
  gas_c = molecular gas

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

initial_values = []
equations = []
r = radio


equations[0] = -(Kh1(r) + Kh2(r)) * gas_H_n(r) - f(r)* gas_H(r) + Wh(r)
equations[1] = Kh1(r) * gas_H_n(r) - D1h(r)
equations[2] = Kh2(r) * gas_H_n(r) - D2h(r)



def gas_H_n(r):
    n = 1.5
    return gas_H(r) ** n