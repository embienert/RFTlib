import numpy as np


def orbit_velocity(mu: float, r: float, a: float) -> float:
    """
    Calculate the current velocity of an object in orbit around a reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param r: [km] Current radius (distance) of the orbiting object from the reference object
    :param a: [km] Semi-major axis of the object's orbit around the reference object
    :return: v [km/s] Current velocity of the object
    """

    return np.sqrt(mu * ((2 / r) - (1 / a)))


def orbit_velocity_circular(mu: float, r: float) -> float:
    """
    Calculate the (constant) velocity of an object in a circular orbit around a reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param r: [km] Radius of the orbit around the reference object
    :return: v [km/s] Velocity of the object on the orbit
    """

    return np.sqrt(mu / r)


def orbit_radius_v_a(mu: float, v: float, a: float) -> float:
    """
    Calculate the current radius (distance) of an object in orbit to its reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param v: [km/s] Current velocity of the orbiting object
    :param a: [km] Semi-major axis of the object's orbit around the reference object
    :return: r [km] Current radius (distance) of the orbiting object from the reference object
    """

    return 1 / (v**2 / (2 * mu) + 1/a)


def orbit_radius_circular_v(mu: float, v: float) -> float:
    """
    Calculate the (constant) radius of an object's circular orbit around a reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param v: [km/s] Velocity of the orbiting object
    :return: r [km] Constant radius (distance) of the orbiting object from the reference object
    """

    return mu / v**2


def orbit_mu_v_r_a(v: float, r: float, a: float) -> float:
    """
    Calculate the standard gravitational parameter of the reference object from an objects orbiting attributes

    :param v: [km/s] Current velocity of the orbiting object
    :param r: [km] Current radius (distance) of the orbiting object from the reference object
    :param a: [km] Semi-major axis of the object's orbit around the reference object
    :return: mu [km³/s²] Standard gravitational parameter of the reference object
    """

    return v**2 / ((2 / r) - (1 / a))


def orbit_mu_circular_v_r(v: float, r: float) -> float:
    """
    Calculate the standard gravitational parameter of the reference object from an objects orbiting attributes,
    assuming the object's orbit is circular (r=a)

    :param v: [km/s] Velocity of the orbiting object
    :param r: [km] Radius of the orbit around the reference object
    :return: mu [km³/s²] Standard gravitational parameter of the reference object
    """

    return v**2 * r


def orbit_period_circular_mu(mu: float, r: float) -> float:
    """
    Calculate the orbit period of an object around a reference object, assuming the object's orbit is circular (r=a)

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param r: [km] Radius of the orbit around the reference object
    :return: T [s] Orbit period of the object
    """

    return 2 * np.pi * np.sqrt(r**3 / mu)


def orbit_period_circular_v(r: float, v: float) -> float:
    """
    Calculate the orbit period of an object around a reference object, assuming the object's orbit is circular (r=a)

    :param r: [km] Radius of the orbit around the reference object
    :param v: [km/s] Velocity of the orbiting object
    :return: T [s] Orbit period of the object
    """

    return 2 * np.pi * r / v


# TODO: orbit period for non-circular orbits

def orbit_radius_T(mu: float, T: float) -> float:
    """
    Calculate the orbit radius of an object around a reference object, assuming the object's orbit is circular (r=a)

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param T: [s] Orbit period of the object
    :return: r [km] Radius of the orbit around the reference object
    """

    return float(np.cbrt(T**2 * mu / (4 * np.pi**2)))


def orbit_velocity_T(r: float, T: float) -> float:
    """
    Calculate the orbit velocity of an object around a reference object, assuming the object's orbit is circular (r=a)

    :param r: [km] Radius of the orbit around the reference object
    :param T: [s] Orbit period of the object
    :return: v [km/s] Velocity of the orbiting object
    """

    return 2 * np.pi * r / T


def average_movement_mu_a(mu: float, a: float) -> float:
    """
    Calculate the average movement ("Mittlere Bewegung") of an object orbiting around a reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param a: Semi-major axis of the object's orbit around the reference object
    :return: n [rad/s] Average movement of the object around the reference object
    """

    return np.sqrt(mu / a**3)


def average_movement_T(T: float) -> float:
    """
    Calculate the average movement ("Mittlere Bewegung") of an object orbiting around a reference object

    :param T: [s] Orbit period of the object
    :return: n [rad/s] Average movement of the object around the reference object
    """

    return 2 * np.pi / T
