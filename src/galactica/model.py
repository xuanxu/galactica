import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from galaxy import Galaxy


g = Galaxy()

time_span = [0, 13.2]
time_eval = np.linspace(time_span[0], time_span[1], 1000)
radius = 1
params = g.parameters(radius)

solution = solve_ivp(g.evolution_step, time_span, g.initial_values(), args=(params,), t_eval=time_eval, method='Radau')
t = solution.t
y = solution.y

sfr_halo, sfr_disk = g.star_formation_rates(y, params)


fig, (plot_main, plot_sf) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

plot_main.set_title("Galaxy evolution")
plot_main.plot(t, y[0], label='Diffuse gas in disk')
plot_main.plot(t, y[1], label='Molecular gas in disk')
plot_main.plot(t, y[2], label='Low mass stars in disk')
plot_main.plot(t, y[3], label='Massive stars in disk')
plot_main.plot(t, y[4], label='Remnants in disk')
plot_main.plot(t, y[5], label='Diffuse gas in halo')
plot_main.plot(t, y[6], label='Low mass stars in halo')
plot_main.plot(t, y[7], label='Massive stars in halo')
plot_main.plot(t, y[8], label='Remnants in halo')
plot_main.set_xlabel("Time (Gyr)")
plot_main.set_ylabel("Mass (Msun)")
plot_main.legend(fontsize=8)
plot_main.grid()

plot_sf.set_title("Star Formation History")
plot_sf.plot(t, sfr_halo, label='SFR halo (Ψ_H)')
plot_sf.plot(t, sfr_disk, label='SFR disk (Ψ_D)')
plot_sf.set_xlabel("Time (Gyr)")
plot_sf.set_ylabel("SFR (Msun / time unit)")
plot_sf.legend()
plot_sf.grid()


plt.tight_layout()
plt.show()
