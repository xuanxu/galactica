import collections
import settings
import numpy as np
import matplotlib.pyplot as plt


class Galaxy:
    def __init__(self, initial_settings=settings.default_settings()):
        self.params = initial_settings

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

    def parameters(self, radius=1):
        params = collections.defaultdict(float)
        params["radius"] = radius
        params["infall_rate"] = 0.1
        params["cloud_formation_factor"] = self.cloud_formation_factor()
        params["star_formation_in_halo"] = self.star_formation_in_halo_factors()  # Kh1 + Kh2
        params["star_formation_cloud_cloud_collisions"] = self.star_formation_cloud_cloud_collisions_factors()  # Ks1 + Ks2 + Ks_restitution_of_diffuse_gas
        params["star_formation_cloud_massive_stars_collisions"] = self.star_formation_cloud_massive_stars_collisions_factors()  # Ka1 + Ka2 + Ka_restitution_of_diffuse_gas

        return params

    def model(self, t, values, params):

        disk_H_gas, disk_H2_gas, disk_low_mass_stars, disk_massive_stars, disk_stars_remnants, \
            halo_H_gas, halo_low_mass_stars, halo_massive_stars, halo_stars_remnants = values

        radius = params["radius"]

        n = 1.5
        disk_H_gas_n = disk_H_gas ** n
        disk_H2_gas_2 = disk_H2_gas ** 2
        halo_H_gas_n = halo_H_gas ** n

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
        Wd = 0
        Wh = 0

        # Death rates for low-mass (<8 Msun) and massive (>8 Msun) stars
        D1d = 0
        D2d = 0
        D1h = 0
        D2h = 0

        infall_rate = params["infall_rate"]

        dt_disk_H_gas = (-Kc * disk_H_gas_n) + (Ka_rest * disk_H2_gas * S2d) + (Ks_rest * disk_H2_gas_2) + (infall_rate * halo_H_gas) + Wd
        dt_disk_H2_gas = (Kc * disk_H_gas_n) - ((Ka1 + Ka2 + Ka_rest) * disk_H2_gas * S2d) - ((Ks1 + Ks2 + Ks_rest) * disk_H2_gas_2)
        dt_disk_low_mass_stars = (Ks1 * disk_H2_gas_2) + (Ka1 * disk_H2_gas * S2d) - D1d
        dt_disk_massive_stars = (Ks2 * disk_H2_gas_2) + (Ka2 * disk_H2_gas * S2d) - D2d
        dt_disk_stars_remnants = D1d + D2d - Wd

        dt_halo_H_gas = -((Kh1 + Kh2) * halo_H_gas_n) - (infall_rate * halo_H_gas) + Wh
        dt_halo_low_mass_stars = (Kh1 * halo_H_gas_n) - D1h
        dt_halo_massive_stars = (Kh2 * halo_H_gas_n) - D2h
        dt_halo_stars_remnants = D1h + D2h - Wh

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
