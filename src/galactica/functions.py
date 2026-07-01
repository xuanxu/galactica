import math
import settings


halo_virial_mass = default_settings()["halo_virial_mass"]
disk_mass = default_settings()["disk_mass"]


def total_mass(radio_kpc):
    return halo_virial_mass + disk_mass * radio_kpc


def disk_mass_at_present_time(radio_kpc):
    return radio_kpc * disk_mass


def delta_total_mass(radio_kpc):
    return total_mass(radio_kpc) - total_mass(radio_kpc-1)


def delta_disk_mass(radio_kpc):
    return disk_mass_at_present_time(radio_kpc) - disk_mass_at_present_time(radio_kpc-1)


def collapse_time(radio_kpc):
    """ infall rate in Gyr """
    tau = -13.2 / math.log(1 - (delta_disk_mass(radio_kpc) / delta_total_mass(radio_kpc)))
    return tau
