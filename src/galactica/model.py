import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from galaxy import Galaxy


g = Galaxy()

time_span = [0, 13.2]
time_eval = np.linspace(time_span[0], time_span[1], 10)
radius = 1

solution = solve_ivp(g.model, time_span, g.initial_values(), args=(g.parameters(radius),), t_eval=time_eval)
