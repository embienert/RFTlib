import numpy as np

from .hohmann import hohmann_circular_velocities_peri, hohmann_circular_velocities_apo, \
    hohmann_velocities_peri, hohmann_velocities_apo


def combined_hohmann_inclination_circular_peri(mu: float, peri_radius: float,
                                               transfer_a: float, delta_i: float) -> float:
    """
    Calculate the velocity required in the peri-center of the origin orbit to transition to the Hohmann orbit with
    changed inclination, assuming both origin and target orbit are circular

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param peri_radius: [km] Radius (distance) of the object in the origin orbit's peri-center
    :param transfer_a: [km] Semi-major axis of the transfer orbit
    :param delta_i: [rad] Inclination change from origin to transfer orbit
    :return: dv12p [km/s] Velocity required for transition to Hohmann orbit with changed inclination
    """

    v1, v2 = hohmann_circular_velocities_peri(mu, peri_radius, transfer_a)
    return np.sqrt(v1**2 + v2**2 - 2 * v1 * v2 * np.cos(delta_i))


def combined_hohmann_inclination_circular_apo(mu: float, transfer_a: float,
                                              apo_radius: float, delta_i: float) -> float:
    """
    Calculate the velocity required in the apo-center of the Hohmann orbit to transition to the target orbit with
    changed inclination, assuming both origin and target orbit are circular

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param transfer_a: [km] Semi-major axis of the transfer orbit
    :param apo_radius: [km] Radius (distance) of the object in the target orbit's apo-center
    :param delta_i: [rad] Inclination change from transfer to target orbit
    :return: dv2a3 [km/s] Velocity required for transition to target orbit with changed inclination
    """

    v1, v2 = hohmann_circular_velocities_apo(mu, transfer_a, apo_radius)
    return np.sqrt(v1**2 + v2**2 - 2 * v1 * v2 * np.cos(delta_i))


def combined_hohmann_inclination_peri(mu: float, origin_a: float, peri_radius: float, transfer_a: float,
                                      delta_i: float) -> float:
    """
    Calculate the velocity required in the peri-center of the origin orbit to transition to the Hohmann orbit with
    changed inclination

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param origin_a: [km] Semi-major axis of the origin orbit
    :param peri_radius: [km] Radius (distance) of the object in the origin orbit's peri-center
    :param transfer_a: [km] Semi-major axis of the transfer orbit
    :param delta_i: [rad] Inclination change from origin to transfer orbit
    :return: dv12p [km/s] Velocity required for transition to Hohmann orbit with changed inclination
    """

    v1, v2 = hohmann_velocities_peri(mu, origin_a, peri_radius, transfer_a)
    return np.sqrt(v1**2 + v2**2 - 2 * v1 * v2 * np.cos(delta_i))


def combined_hohmann_inclination_apo(mu: float, transfer_a: float,
                                     target_a: float, apo_radius: float, delta_i: float) -> float:
    """
    Calculate the velocity required in the apo-center of the Hohmann orbit to transition to the target orbit with
    changed inclination

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param transfer_a: [km] Semi-major axis of the transfer orbit
    :param target_a: [km] Semi-major axis of the target orbit
    :param apo_radius: [km] Radius (distance) of the object in the target orbit's apo-center
    :param delta_i: [rad] Inclination change from transfer to target orbit
    :return: dv2a3 [km/s] Velocity required for transition to target orbit with changed inclination
    """

    v1, v2 = hohmann_velocities_apo(mu, transfer_a, target_a, apo_radius)
    return np.sqrt(v1**2 + v2**2 - 2 * v1 * v2 * np.cos(delta_i))
