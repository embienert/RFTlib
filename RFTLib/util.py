import numpy as np


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
