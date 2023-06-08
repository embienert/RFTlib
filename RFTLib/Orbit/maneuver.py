import numpy as np

from .orbit import orbit_velocity, orbit_velocity_circular


# ------------------------------------ #
# Inclination change w/o radius change #
# ------------------------------------ #

def inclination_change(delta_i: float, v: float):
    """
    Calculate required velocity for an inclination maneuver of the object's current orbit

    :param delta_i: [rad] Desired inclination maneuver
    :param v: [km/s] Current velocity of the object
    :return: v [km/s] Required velocity for the inclination maneuver
    """

    return 2 * v * np.sin(delta_i / 2)


def inclination_inclination_change_dv(dv: float, v: float):
    """
    Calculate the inclination change

    :param dv: [km/s] Velocity change of the inclination maneuver
    :param v: [km/s] Current Velocity of the object
    :return: delta_i [] Inclination change
    """

    return 2 * np.arcsin(dv / (2 * v))


# ------------------------------------ #
# Hohmann transfer to different radius #
# ------------------------------------ #

def hohmann_transfer_a(origin_radius: float, target_radius: float) -> float:
    """
    Calculate the semi-major axis of the hohmann transfer orbit

    :param origin_radius: [km] Radius of the peri-center of the origin orbit
    :param target_radius: [km] Radius of the apo-center of the target orbit
    :return: a [km] Semi-major axis of the hohmann transfer orbit
    """

    return (origin_radius + target_radius) / 2


def hohmann_circular_velocities_peri(mu: float, peri_radius: float, transfer_a: float) -> (float, float):
    """
    Calculate the velocity of an object during a Hohmann transfer in the peri-center of its origin and transfer orbit
    around a reference object, assuming both origin and target orbit are circular

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param peri_radius: [km] Radius (distance) of the object in the origin orbit's peri-center
    :param transfer_a: [km] Semi-major axis of the transfer orbit
    :return: (v1, v2p) ([km/s], [km/s]) Velocity in the origin and transfer orbit's peri-center
    """

    origin_peri_velocity = orbit_velocity_circular(mu, peri_radius)
    transfer_peri_velocity = orbit_velocity(mu, peri_radius, transfer_a)

    return origin_peri_velocity, transfer_peri_velocity


def hohmann_circular_velocities_apo(mu: float, transfer_a: float, apo_radius: float) -> (float, float):
    """
    Calculate the velocity of an object during a Hohmann transfer in the apo-center of its transfer and target orbit
    around a reference object, assuming both origin and target orbit are circular

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param transfer_a: [km] Semi-major axis of the transfer orbit
    :param apo_radius: [km] Radius (distance) of the object in the target orbit's apo-center
    :return: (v2a, v3) ([km/s], [km/s]) Velocity in the transfer and target orbit's apo-center
    """

    transfer_apo_velocity = orbit_velocity(mu, apo_radius, transfer_a)
    target_apo_velocity = orbit_velocity_circular(mu, apo_radius)

    return transfer_apo_velocity, target_apo_velocity


def hohmann_circular_velocities(mu: float, peri_radius: float, apo_radius: float) -> (float, float, float, float):
    """
    Calculate the velocities of origin, transfer and target orbit in both peri- and apo-center, assuming both origin
    and target orbit are circular

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param peri_radius: [km] Radius (distance) of the object in the origin orbit's peri-center
    :param apo_radius: [km] Radius (distance) of the object in the target orbit's apo-center
    :return: (v1, v2a, v2p, v3) ([km/s], [km/s], [km/s], [km/s]) Velocity in the origin and transfer orbit's
        peri-center, as well as the in the transfer and target orbit's apo-center
    """

    transfer_a = hohmann_transfer_a(peri_radius, apo_radius)

    return *hohmann_circular_velocities_peri(mu, peri_radius, transfer_a), \
        *hohmann_circular_velocities_apo(mu, transfer_a, apo_radius)


def hohmann_velocities_peri(mu: float, origin_a: float, peri_radius: float, transfer_a: float) -> (float, float):
    """
    Calculate the velocity of an object during a Hohmann transfer in the peri-center of its origin and transfer orbit
    around a reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param origin_a: [km] Semi-major axis of the origin orbit
    :param peri_radius: [km] Radius (distance) of the object in the origin orbit's peri-center
    :param transfer_a: [km] Semi-major axis of the transfer orbit
    :return: (v1, v2p) ([km/s], [km/s]) Velocity in the origin and transfer orbit's peri-center
    """

    origin_peri_velocity = orbit_velocity(mu, peri_radius, origin_a)
    transfer_peri_velocity = orbit_velocity(mu, peri_radius, transfer_a)

    return origin_peri_velocity, transfer_peri_velocity


def hohmann_velocities_apo(mu: float, transfer_a: float, target_a: float, apo_radius: float) -> (float, float):
    """
    Calculate the velocity of an object during a Hohmann transfer in the apo-center of its transfer and target orbit
    around a reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param transfer_a: [km] Semi-major axis of the transfer orbit
    :param target_a: [km] Semi-major axis of the target orbit
    :param apo_radius: [km] Radius (distance) of the object in the target orbit's apo-center
    :return: (v2a, v3) ([km/s], [km/s]) Velocity in the transfer and target orbit's apo-center
    """

    transfer_apo_velocity = orbit_velocity(mu, apo_radius, transfer_a)
    target_apo_velocity = orbit_velocity(mu, apo_radius, target_a)

    return transfer_apo_velocity, target_apo_velocity


def hohmann_velocities(mu: float, origin_a: float, peri_radius: float,
                       target_a: float, apo_radius: float) -> (float, float, float, float):
    """
    Calculate the velocities of origin, transfer and target orbit in both peri- and apo-center

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param origin_a: [km] Semi-major axis of the origin orbit
    :param peri_radius: [km] Radius (distance) of the object in the origin orbit's peri-center
    :param target_a: [km] Semi-major axis of the target orbit
    :param apo_radius: [km] Radius (distance) of the object in the target orbit's apo-center
    :return: (v1, v2a, v2p, v3) ([km/s], [km/s], [km/s], [km/s]) Velocity in the origin and transfer orbit's
        peri-center, as well as the in the transfer and target orbit's apo-center
    """

    transfer_a = hohmann_transfer_a(peri_radius, apo_radius)

    return *hohmann_velocities_peri(mu, origin_a, peri_radius, transfer_a), \
        *hohmann_velocities_apo(mu, transfer_a, target_a, apo_radius)


# --------------------------------------------------------------------------------------------------- #
# Hohmann transfer to/from transfer orbit with simultaneous inclination change in peri- or apo-center #
# --------------------------------------------------------------------------------------------------- #

def combined_hohmann_inclination_circular_peri(mu: float, peri_radius: float,
                                               transfer_a: float, delta_i: float) -> float:
    """
    Calculate the velocity required in the peri-center of the origin orbit to transition to the Hohmann orbit with
    changed inclination, assuming both origin and target orbit are circular

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param peri_radius: [km] Radius (distance) of the object in the origin orbit's peri-center
    :param transfer_a: [km] Semi-major axis of the transfer orbit
    :param delta_i: [rad] Inclination maneuver from origin to transfer orbit
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
    :param delta_i: [rad] Inclination maneuver from transfer to target orbit
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
    :param delta_i: [rad] Inclination maneuver from origin to transfer orbit
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
    :param delta_i: [rad] Inclination maneuver from transfer to target orbit
    :return: dv2a3 [km/s] Velocity required for transition to target orbit with changed inclination
    """

    v1, v2 = hohmann_velocities_apo(mu, transfer_a, target_a, apo_radius)
    return np.sqrt(v1**2 + v2**2 - 2 * v1 * v2 * np.cos(delta_i))
