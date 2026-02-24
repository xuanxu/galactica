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
