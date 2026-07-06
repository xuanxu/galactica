import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from galactica.galaxy import Galaxy


def precompute_death_rates(g, t_grid, psi_halo_interp, psi_disk_interp):
    D1h = np.array([g.death_rate(t, psi_halo_interp, 0.15, 4.0) for t in t_grid])
    D2h = np.array([g.death_rate(t, psi_halo_interp, 4.0, 100.0) for t in t_grid])
    D1d = np.array([g.death_rate(t, psi_disk_interp, 0.15, 4.0) for t in t_grid])
    D2d = np.array([g.death_rate(t, psi_disk_interp, 4.0, 100.0) for t in t_grid])
    kw = dict(bounds_error=False, fill_value=0.0)
    return (interp1d(t_grid, D1h, **kw), interp1d(t_grid, D2h, **kw),
            interp1d(t_grid, D1d, **kw), interp1d(t_grid, D2d, **kw))


def run_model(g, params, t_span, t_eval, death_rate_interps=None):
  solution = solve_ivp(g.evolution_step, time_span, g.initial_values(), args=(params,), t_eval=time_eval, method='Radau')
  sfr_halo, sfr_disk = g.star_formation_rates(solution.y, params)
  return solution, sfr_halo, sfr_disk


g = Galaxy()

time_span = [0, 13.2]
time_eval = np.linspace(time_span[0], time_span[1], 1000)
radius = 8
params = g.parameters(radius)


solution, sfr_halo, sfr_disk = None, None, None
death_rates_interps = None

for pass_num in range(1, 4):
    print(f"Pass {pass_num}/3...")
    if solution is not None:
        psi_halo_interp = interp1d(solution.t, sfr_halo, bounds_error=False, fill_value=0.0)
        psi_disk_interp = interp1d(solution.t, sfr_disk, bounds_error=False, fill_value=0.0)
        death_rates_interps = precompute_death_rates(g, time_eval, psi_halo_interp, psi_disk_interp)
    solution, sfr_halo, sfr_disk = run_model(g, params, time_span, time_eval, death_rates_interps)

t = solution.t
y = solution.y

fig, (plot_disk, plot_halo, plot_sf) = plt.subplots(3, 1, figsize=(10, 11), sharex=True)

plot_disk.set_title("Disk evolution")
plot_disk.plot(t, y[0], label='Diffuse gas')
plot_disk.plot(t, y[1], label='Molecular gas')
plot_disk.plot(t, y[2], label='Low mass stars')
plot_disk.plot(t, y[3], label='Massive stars')
plot_disk.plot(t, y[4], label='Remnants')
plot_disk.set_xlabel("Time (Gyr)")
plot_disk.set_ylabel("Mass (Msun)")
plot_disk.legend(fontsize=8)
plot_disk.grid()

plot_halo.set_title("Halo evolution")
plot_halo.plot(t, y[5], label='Diffuse gas')
plot_halo.plot(t, y[6], label='Low mass stars')
plot_halo.plot(t, y[7], label='Massive stars')
plot_halo.plot(t, y[8], label='Remnants')
plot_halo.set_xlabel("Time (Gyr)")
plot_halo.set_ylabel("Mass (Msun)")
plot_halo.legend(fontsize=8)
plot_halo.grid()

plot_sf.set_title("Star Formation History")
plot_sf.plot(t, sfr_halo, label='SFR halo (Ψ_H)')
plot_sf.plot(t, sfr_disk, label='SFR disk (Ψ_D)')
plot_sf.set_xlabel("Time (Gyr)")
plot_sf.set_ylabel("SFR (Msun / time unit)")
plot_sf.legend()
plot_sf.grid()


plt.tight_layout()
plt.show()
