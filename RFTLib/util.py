import numpy as np


# --------- #
# Constants #
# --------- #

G = 6.6743e-20  # km³/(kg*s²)

# EARTH Constants
EARTH_MU = 3.986004415e5  # km³/s²
EARTH_MASS = 5.9722e24  # kg
EARTH_G = 9.81e-3  # km/s²
EARTH_SURFACE_ESCAPE_VELOCITY = 11.183  # km/s
EARTH_SEMIMAJOR_AXIS = 149.6e6  # km
EARTH_RADIUS = 6371  # km
EARTH_DAY = 86400  # s
EARTH_STARDAY = 86164.0916  # s
EARTH_STARYEAR = 365.256_360_42 * (60 * 60 * 24)
EARTH_ROTATION_VELOCITY = 2 * np.pi * EARTH_RADIUS / EARTH_DAY  # km/s
EARTH_ROTATION_PERIOD = 23.9 * 60 * 60  # s
EARTH_J2 = 0.001_082_626_68  # no unit

# MOON Constants
MOON_MU = 4.9048695e3  # km³/s²
MOON_MASS = 7.342e22  # kg
MOON_G = 1.6220e-3  # km/s²
MOON_ESCAPE_VELOCITY = 2.3762  # km/s
MOON_SEMIMAJOR_AXIS = 3.84399e5  # km
MOON_RADIUS = 1738.1  # km
MOON_ROTATION_VELOCITY = 4.627e-3  # km/s

# SUN Constants
SUN_MU = 1.327_124_400_41e11  # km³/s²
SUN_MASS = 1.9884e30  # kg
SUN_G = 274e-3  # km/s²
SUN_ESCAPE_VELOCITY = 617.6  # km/s


def deg2rad(deg: float) -> float:
    """
    Converts degrees to radians

    :param deg: [deg] Value in degrees
    :return: [rad] Value in radians
    """

    return (deg % 360) * np.pi / 180


def rad2deg(rad: float) -> float:
    """
    Converts radians to degrees

    :param rad: [rad] Value in radians
    :return: [deg] Value in degrees
    """

    return (rad % (2 * np.pi)) * 180 / np.pi


def delta(val1: float, val2: float) -> float:
    """
    Calculates the difference between the two values
    """

    return np.abs(val1 - val2)


def s_to_d(val_s: float) -> float:
    """
    Converts seconds to days
    """

    return val_s * (60 * 60 * 24)


def d_to_s(val_d: float) -> float:
    """
    Converts days to seconds
    """

    return val_d / (60 * 60 * 24)
