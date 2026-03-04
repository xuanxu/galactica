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
    defaults["G"] = 0.44985  # Gravitational constant in Kpc^3/(10^9Msun * 10^7yrs)

    return defaults


def default_initial_values():
    # Disk
    disk_H_gas = 1e11
    disk_H2_gas = 0.0
    disk_low_mass_stars = 0.0
    disk_massive_stars = 0.0
    disk_stars_remnants = 0.0
    disk_death_low_mass_stars = 0.0
    disk_death_massive_stars = 0.0
    # Halo
    halo_H_gas = 1e12
    halo_H2_gas = 0.0
    halo_low_mass_stars = 0.0
    halo_massive_stars = 0.0
    halo_stars_remnants = 0.0

    return [disk_H_gas, disk_H2_gas, disk_low_mass_stars, disk_massive_stars, disk_stars_remnants,
            halo_H_gas, halo_low_mass_stars, halo_massive_stars, halo_stars_remnants]
