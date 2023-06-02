import numpy as np

from orbit import average_movement_mu_a, average_movement_T


# TODO: Operations regarding increasing/decreasing knot

def knot_drift_J2_mu(mu: float, r: float, a: float, e: float, i: float, J2: float) -> float:
    """
    Calculate knot drift of an object's orbit around a reference object using J2 model

    :param mu: [km³/s²] Standard gravitational parameter of the reference system
    :param r: [km] Radius of the reference system
    :param a: [km] Semi-major axis of the object's orbit around the reference object
    :param e: [] Numeric eccentricity of the object's orbit around the reference object
    :param i: [rad] Inclination of the object's orbit around the reference object
    :param J2: [] J2 coefficient describing the "oblateness" of the reference object
    :return: Omega(t) [rad/s] Knot drift of the object's orbit around the reference object
    """

    n = average_movement_mu_a(mu, a)

    return -3/2 * n * J2 * (r / a)**2 * np.cos(i) / (1 - e**2)**2


def knot_drift_J2_T(r: float, T: float, a: float, e: float, i: float, J2: float) -> float:
    """
    Calculate knot drift of an object's orbit around a reference object using J2 model

    :param r: [km] Radius of the reference system
    :param T: [s] Orbit period of the object
    :param a: [km] Semi-major axis of the object's orbit around the reference object
    :param e: [] Numeric eccentricity of the object's orbit around the reference object
    :param i: [rad] Inclination of the object's orbit around the reference object
    :param J2: [] J2 coefficient describing the "oblateness" of the reference object
    :return: Omega(t) [rad/s] Knot drift of the object's orbit around the reference object
    """

    n = average_movement_T(T)

    return -3 / 2 * n * J2 * (r / a) ** 2 * np.cos(i) / (1 - e ** 2) ** 2


def inclination_knot_drift_J2_mu(mu: float, r: float, a: float, e: float, drift: float, J2: float) -> float:
    """
    Calculate the required inclination of an object's orbit around a reference object to achieve a desired knot drift
    using the J2 model

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param r: [km] Radius of the reference system
    :param a: [km] Semi-major axis of the object's orbit around the reference object
    :param e: [] Numeric eccentricity of the object's orbit around the reference object
    :param drift: [rad/s] Knot drift of the object's orbit around the reference object
    :param J2: [] J2 coefficient describing the "oblateness" of the reference object
    :return: i [rad] Required inclination of the object's orbit around the reference object
    """

    n = average_movement_mu_a(mu, a)

    return np.arccos(-2 * drift * a ** 2 * (1 - e ** 2) ** 2 / (3 * n * J2 * r ** 2))


def inclination_knot_drift_J2_T(r: float, T: float, a: float, e: float, drift: float, J2: float) -> float:
    """
    Calculate the required inclination of an object's orbit around a reference object to achieve a desired knot drift
    using the J2 model

    :param r: [km] Radius of the reference system
    :param T: [s] Orbit period of the object
    :param a: [km] Semi-major axis of the object's orbit around the reference object
    :param e: [] Numeric eccentricity of the object's orbit around the reference object
    :param drift: [rad/s] Knot drift of the object's orbit around the reference object
    :param J2: [] J2 coefficient describing the "oblateness" of the reference object
    :return: i [rad] Required inclination of the object's orbit around the reference object
    """

    n = average_movement_T(T)

    return np.arccos(-2 * drift * a ** 2 * (1 - e ** 2) ** 2 / (3 * n * J2 * r ** 2))
