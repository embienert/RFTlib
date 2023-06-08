import numpy as np


# --------- #
# Constants #
# --------- #

G = 6.6743e-20  # km^3/(kg*s^2)

# EARTH Constants
EARTH_MU = 3.986004418e5  # km^3/s^2
EARTH_MASS = 5.9722e24  # kg
EARTH_G = 9.806e-3  # km/s^2
EARTH_ESCAPE_VELOCITY = 11.183  # km/s
EARTH_SEMIMAJOR_AXIS = 149.6e6  # km
EARTH_RADIUS = 6371  # km

# SUN Constants
SUN_MU = 1.327_124_400_41e11  # km^3/s^2
SUN_MASS = 1.9884e30  # kg
SUN_G = 274e-3  # km/s^2
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

    return np.abs(val1, val2)
