import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def disk_gas(t, y, k):


  Kc = 1 # star_formation_factor_cloud
  gas_D = 1e11
  n = 1.5
  gas_D_n = gas_D ** n

  return(-Kc * gas_D_n) + (Ka_rest * c * S2d) + (Ks_rest * c2) + (f * gas_H) + Wd