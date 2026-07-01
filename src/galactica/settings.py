"""
Default settings and initial parameters

"""

import collections


def default_settings():
    defaults = collections.defaultdict(float)
    defaults["region_width_kpc"] = 1.0  # width of a default region in kiloparsecs
    defaults["region_galactocentric_radio_kpc"] = 8.0  # distance of the region to the galactic center in kiloparsecs
    defaults["halo_radio_kpc"] = 260.996  # radio of the halo in kiloparsecs
    defaults["disk_height_kpc"] = 0.2  # height of the galactic disk in kiloparsecs
    defaults["G"] = 0.44985   # Gravitational constant in Kpc^3/(10^9Msun * 10^7yrs)
    defaults["boundary_mass"] = 4.0  # Msun, boundary between low-mass and massive star groups
    defaults["halo_virial_mass"] = 1e12  # [5e10, 1e13] Msol
    defaults["disk_mass"] = 1e11  # [1.25e8, 5.3e11] Msol
    defaults["max_radio"] = 16  # kpc
    defaults["metallicity"] = 0.02  # solar metallicity
    defaults["imf_min_mass"] = 0.1  # Msun
    defaults["imf_max_mass"] = 100  # Msun

    return defaults


def default_initial_values():
    # Disk
    disk_H_gas = default_settings()["disk_mass"]
    disk_H2_gas = 0.0
    disk_low_mass_stars = 0.0
    disk_massive_stars = 0.0
    disk_stars_remnants = 0.0
    disk_death_low_mass_stars = 0.0
    disk_death_massive_stars = 0.0
    # Halo
    halo_H_gas = default_settings()["halo_virial_mass"]
    halo_H2_gas = 0.0
    halo_low_mass_stars = 0.0
    halo_massive_stars = 0.0
    halo_stars_remnants = 0.0

    return [disk_H_gas, disk_H2_gas, disk_low_mass_stars, disk_massive_stars, disk_stars_remnants,
            halo_H_gas, halo_low_mass_stars, halo_massive_stars, halo_stars_remnants]
