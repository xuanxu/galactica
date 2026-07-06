import collections
import galactica.settings as settings
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.optimize import brentq
from scipy.integrate import quad, IntegrationWarning
import starmatrix
from starmatrix import elements as sm_elements
from starmatrix import abundances as sm_abundances
from starmatrix import imfs as sm_imfs
from starmatrix import functions as sm_functions
from starmatrix import settings as sm_settings


class Galaxy:
    def __init__(self, initial_settings=settings.default_settings()):
        self.params = initial_settings
        self.imf = sm_imfs.select_imf("kroupa2002", {})
        self.expelled = sm_elements.Expelled(starmatrix.settings.default["expelled_elements_filename"])
        self.abundances = sm_abundances.select_abundances("as09", self.params["metallicity"])
        self.return_fraction = sm_functions.return_fraction(self.params["imf_min_mass"], self.params["imf_max_mass"], self.expelled, self.imf)

    def initial_values(self):
        return settings.default_initial_values()

    def cloud_formation_factor(self):
        return 1

    def star_formation_cloud_massive_stars_collisions_factors(self):
        low_stars_factor = 1.0
        massive_stars_factor = 1.0
        restitution_to_diffuse_gas = 0.1
        return [low_stars_factor, massive_stars_factor, restitution_to_diffuse_gas]

    def star_formation_cloud_cloud_collisions_factors(self):
        low_stars_factor = 1.0
        massive_stars_factor = 1.0
        restitution_to_diffuse_gas = 0.1
        return [low_stars_factor, massive_stars_factor, restitution_to_diffuse_gas]

    def star_formation_in_halo_factors(self):
        efficiency = 2.173  # epsilon_h computed for a best value K_h = 9e-3 able to reproduce SFR and abundances of MWG halo
        factor = efficiency * (self.params["G"] / self.volume_halo())**0.5

        low_stars_factor = factor*0.5
        massive_stars_factor = factor*0.5
        return [low_stars_factor, massive_stars_factor]

    def volume_halo(self, region_shape='ring'):
        h = np.sqrt((self.params["halo_radio_kpc"] ** 2) - (self.params["region_galactocentric_radio_kpc"] ** 2))
        if region_shape == 'square':
            square_area = self.params["region_width_kpc"] ** 2
            return square_area * 2 * h
        elif region_shape == 'ring':
            half_ring_width = 0.5 * self.params["region_width_kpc"]
            ring_area = np.pi * (
                (self.params["region_galactocentric_radio_kpc"] + half_ring_width) ** 2 -
                (self.params["region_galactocentric_radio_kpc"] - half_ring_width) ** 2)
            return ring_area * 2 * h
        else:
            raise Exception("Wrong region shape for halo. Allowed options: [square, ring]")

    def volume_disk(self, region_shape='ring'):
        h = self.params["disk_height_kpc"]
        if region_shape == 'square':
            square_area = self.params["region_width_kpc"] ** 2
            return square_area * h
        elif region_shape == 'ring':
            half_ring_width = 0.5 * self.params["region_width_kpc"]
            ring_area = np.pi * (
                (self.params["region_galactocentric_radio_kpc"] + half_ring_width) ** 2 -
                (self.params["region_galactocentric_radio_kpc"] - half_ring_width) ** 2)
            return ring_area * h
        else:
            raise Exception("Wrong region shape for disk. Allowed options: [square, ring]")

    def star_formation_rates(self, values, params):
        """Compute SFR for halo and disk

        ΨH = (Kh1 + Kh2) * gH^n
        ΨD = (Ks1 + Ks2) * c^2 + (Ka1 + Ka2) * c * s2D
        """
        values = np.asarray(values)
        disk_H2_gas = values[1]
        disk_massive_stars = values[3]
        halo_H_gas = values[5]

        Kh1, Kh2 = params["star_formation_in_halo"]
        Ks1, Ks2 = params["star_formation_cloud_cloud_collisions"][:2]
        Ka1, Ka2 = params["star_formation_cloud_massive_stars_collisions"][:2]

        psi_halo = (Kh1 + Kh2) * halo_H_gas ** 1.5
        psi_disk = (Ks1 + Ks2) * disk_H2_gas ** 2 + (Ka1 + Ka2) * disk_H2_gas * disk_massive_stars

        return psi_halo, psi_disk

    def stellar_lifetime(self, m, z=None):
        """Main-sequence lifetime in Gyr via Starmatrix."""
        if z is None:
            z = self.params["metallicity"]

        return sm_functions.stellar_lifetime(m, z)

    def parameters(self, radius=1):
        params = collections.defaultdict(float)
        params["radius"] = radius
        params["infall_rate"] = 0.1
        params["cloud_formation_factor"] = self.cloud_formation_factor()
        params["star_formation_in_halo"] = self.star_formation_in_halo_factors()  # Kh1 + Kh2
        params["star_formation_cloud_cloud_collisions"] = self.star_formation_cloud_cloud_collisions_factors()  # Ks1 + Ks2 + Ks_restitution_of_diffuse_gas
        params["star_formation_cloud_massive_stars_collisions"] = self.star_formation_cloud_massive_stars_collisions_factors()  # Ka1 + Ka2 + Ka_restitution_of_diffuse_gas

        return params

    def death_rate(self, t, psi_interp, m_low, m_high):
        """Death rate integral at time t

        D(t) = ∫[m_crit, m_high] ψ(t − τ(m)) · m · φ(m) dm

        """
        if t <= self.stellar_lifetime(m_high):
            return 0.0

        if t >= self.stellar_lifetime(m_low):
            m_crit = m_low
        else:
            m_crit = brentq(lambda m: self.stellar_lifetime(m) - t, m_low, m_high, xtol=1e-4)

        def integrand(m):
            tau = self.stellar_lifetime(m)
            return float(psi_interp(t - tau)) * m * self.imf.for_mass(m)

        # Hint quad about Kroupa2002 IMF slope changes at 0.08, 0.5, 1.0 Msun
        pts = [k for k in [0.08, 0.5, 1.0] if m_crit < k < m_high]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            result, _ = quad(integrand, m_crit, m_high, limit=200, points=pts)
        return result

    def evolution_step(self, t, values, params):

        disk_H_gas, disk_H2_gas, disk_low_mass_stars, disk_massive_stars, disk_stars_remnants, \
            halo_H_gas, halo_low_mass_stars, halo_massive_stars, halo_stars_remnants = values

        radius = params["radius"]

        n = 1.5
        disk_H_gas_n = max(0.0, disk_H_gas) ** n
        disk_H2_gas_2 = max(0.0, disk_H2_gas) ** 2
        halo_H_gas_n = max(0.0, halo_H_gas) ** n

        # Star Formation factors
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

        S2d = 1.0  # s_massive_disk dummy value

        # Death rates for low-mass (<8 Msun) and massive (>8 Msun) stars
        disk_death_rate_low_mass = 0
        disk_death_rate_massive = 0
        halo_death_rate_low_mass = 0
        halo_death_rate_massive = 0

        D1h_interp = params.get("D1h_interp")
        if D1h_interp is not None:
            disk_death_rate_low_mass = float(params["D1d_interp"](t))
            disk_death_rate_massive = float(params["D2d_interp"](t))
            halo_death_rate_low_mass = float(params["D1h_interp"](t))
            halo_death_rate_massive = float(params["D2h_interp"](t))

        Wd = self.return_fraction * (disk_death_rate_low_mass + disk_death_rate_massive)
        Wh = self.return_fraction * (halo_death_rate_low_mass + halo_death_rate_massive)

        infall_rate = params["infall_rate"]

        dt_disk_H_gas = (-Kc * disk_H_gas_n) + (Ka_rest * disk_H2_gas * S2d) + (Ks_rest * disk_H2_gas_2) + (infall_rate * halo_H_gas) + Wd
        dt_disk_H2_gas = (Kc * disk_H_gas_n) - ((Ka1 + Ka2 + Ka_rest) * disk_H2_gas * S2d) - ((Ks1 + Ks2 + Ks_rest) * disk_H2_gas_2)
        dt_disk_low_mass_stars = (Ks1 * disk_H2_gas_2) + (Ka1 * disk_H2_gas * S2d) - disk_death_rate_low_mass
        dt_disk_massive_stars = (Ks2 * disk_H2_gas_2) + (Ka2 * disk_H2_gas * S2d) - disk_death_rate_massive
        dt_disk_stars_remnants = disk_death_rate_low_mass + disk_death_rate_massive - Wd

        dt_halo_H_gas = -((Kh1 + Kh2) * halo_H_gas_n) - (infall_rate * halo_H_gas) + Wh
        dt_halo_low_mass_stars = (Kh1 * halo_H_gas_n) - halo_death_rate_low_mass
        dt_halo_massive_stars = (Kh2 * halo_H_gas_n) - halo_death_rate_massive
        dt_halo_stars_remnants = halo_death_rate_low_mass + halo_death_rate_massive - Wh

        dtdy = [dt_disk_H_gas,
                dt_disk_H2_gas,
                dt_disk_low_mass_stars,
                dt_disk_massive_stars,
                dt_disk_stars_remnants,
                dt_halo_H_gas,
                dt_halo_low_mass_stars,
                dt_halo_massive_stars,
                dt_halo_stars_remnants]

        return dtdy
