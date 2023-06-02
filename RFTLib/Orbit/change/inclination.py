import numpy as np


def inclination_change(delta_i: float, v: float):
    """
    Calculate required velocity for an inclination change of the object's current orbit

    :param delta_i: [rad] Desired inclination change
    :param v: [km/s] Current velocity of the object
    :return: v [km/s] Required velocity for the inclination change
    """

    return 2 * v * np.sin(delta_i / 2)
