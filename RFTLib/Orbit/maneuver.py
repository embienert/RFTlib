import numpy as np

from .orbit import orbit_velocity, orbit_velocity_circular, vesc_r, v_vesc
from ..util import delta


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
    return np.sqrt(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * np.cos(delta_i))


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
    return np.sqrt(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * np.cos(delta_i))


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
    return np.sqrt(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * np.cos(delta_i))


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
    return np.sqrt(v1 ** 2 + v2 ** 2 - 2 * v1 * v2 * np.cos(delta_i))


def ascending_knot_inclination(v: float, i_origin: float, i_target: float, delta_Omega: float) -> (float, float):
    """
    Calculate the required velocity and length to change inclination and ascending knot of the orbit

    :param v: [km/s] Current velocity on the orbit
    :param i_origin: [rad] Inclination of the origin orbit
    :param i_target: [rad] Inclination of the target orbit
    :param delta_Omega: [rad] Difference between origin and target ascending knot
    :return: (dv, lambda) ([km/s], rad) Required velocity and length to change inclination and ascending knot
        of the orbit
    """

    alpha = np.arccos(np.cos(i_origin) * np.cos(i_target) + np.sin(i_origin) * np.sin(i_target) * np.cos(delta_Omega))
    lam = np.sin(i_target) * np.sin(delta_Omega) / np.sin(alpha)

    dv = 2 * v * np.sin(alpha / 2)

    return dv, lam


def interplanetary_hohmann_circular(origin_mu: float, target_mu: float, reference_mu: float,
                                    reference_origin_radius: float, reference_target_radius: float,
                                    origin_radius: float, target_radius: float) -> (float, float):
    """
    Calculate the required velocity differences to transfer an object from its origin orbit in the origin system
    to the target orbit in the target system, assuming all the discussed orbits are circular

    :param origin_mu: [km³/s²] Standard gravitational parameter of the origin system
    :param target_mu: [km³/s²] Standard gravitational parameter of the target system
    :param reference_mu: [km³/s²] Standard gravitational parameter of the reference system
    :param reference_origin_radius: [km] Orbit radius of the origin system around the reference system
    :param reference_target_radius: [km] Orbit radius of the target system around the reference system
    :param origin_radius: [km] Orbit radius of the object around the origin system
    :param target_radius: [km] Orbit radius of the object around the target system
    :return: (dv12p, dv4p5) ([km/s], [km/s]) Required velocity differences in peri-center of the origin and target orbit
    """

    # Orbit velocities of the origin and target system around the reference system
    orbit_velocity_origin_system = orbit_velocity_circular(reference_mu, reference_origin_radius)
    orbit_velocity_target_system = orbit_velocity_circular(reference_mu, reference_target_radius)

    # Semi-major axis of the transfer orbit
    transfer_a = (reference_origin_radius + reference_target_radius) / 2

    # Required velocity in the peri- and apo-center of the transfer orbit
    v3p = orbit_velocity(reference_mu, reference_origin_radius, transfer_a)
    v3a = orbit_velocity(reference_mu, reference_target_radius, transfer_a)

    # Calculate the hyperbolic excess velocities at SOI edge of origin and target system
    hyperbolic_v_inf_origin = delta(v3p, orbit_velocity_origin_system)
    hyperbolic_v_inf_target = delta(v3a, orbit_velocity_target_system)

    # Calculate escape velocities of origin and target orbit within their system
    origin_v_esc = vesc_r(origin_mu, origin_radius)
    target_v_esc = vesc_r(target_mu, target_radius)

    # Velocity in the peri-center of the parabolic orbits within their system
    v2p = v_vesc(origin_v_esc, hyperbolic_v_inf_origin)
    v4p = v_vesc(target_v_esc, hyperbolic_v_inf_target)

    # Calculate the velocities of the origin and target orbits within their system
    v1 = orbit_velocity_circular(origin_mu, origin_radius)
    v5 = orbit_velocity_circular(target_mu, target_radius)

    # Calculate the required velocity differences in the peri-center of origin and target orbit within their systems
    dv12p = delta(v1, v2p)
    dv4p5 = delta(v4p, v5)

    return dv12p, dv4p5


def interplanetary_hohmann(origin_mu: float, target_mu: float, reference_mu: float,
                           reference_origin_radius: float, reference_origin_a: float,
                           reference_target_radius: float, reference_target_a: float,
                           origin_radius: float, origin_a: float,
                           target_radius: float, target_a: float) -> (float, float):
    """
    Calculate the required velocity differences to transfer an objects from its origin orbit in the origin system
    to the target orbit in the target system

    :param origin_mu: [km³/s²] Standard gravitational parameter of the origin system
    :param target_mu: [km³/s²] Standard gravitational parameter of the target system
    :param reference_mu: [km³/s²] Standard gravitational parameter of the reference system
    :param reference_origin_radius: [km] Current distance of the origin system to the reference system
    :param reference_origin_a: [km] Semi-major axis of the origin system's orbit around the reference system
    :param reference_target_radius: [km] Current distance of the target system to the reference system
    :param reference_target_a: [km] Semi-major axis of the target system's orbit around the reference system
    :param origin_radius: [km] Current distance of the object to the origin system
    :param origin_a: [km] Semi-major axis of the object's orbit around the origin system
    :param target_radius: [km] Current distance of the object to the target system
    :param target_a: [km] Semi-major axis of the object's orbit around the target system
    :return: (dv12p, dv4p5) ([km/s], [km/s]) Required velocity differences in peri-center of the origin and target orbit
    """

    # Orbit velocities of the origin and target system around the reference system
    orbit_velocity_origin_system = orbit_velocity(reference_mu, reference_origin_radius, reference_origin_a)
    orbit_velocity_target_system = orbit_velocity(reference_mu, reference_target_radius, reference_target_a)

    # Semi-major axis of the transfer orbit
    transfer_a = (reference_origin_a + reference_target_a) / 2

    # Required velocity in the peri- and apo-center of the transfer orbit
    v3p = orbit_velocity(reference_mu, reference_origin_radius, transfer_a)
    v3a = orbit_velocity(reference_mu, reference_target_radius, transfer_a)

    # Calculate the hyperbolic excess velocities at SOI edge of origin and target system
    hyperbolic_v_inf_origin = delta(v3p, orbit_velocity_origin_system)
    hyperbolic_v_inf_target = delta(v3a, orbit_velocity_target_system)

    # Calculate escape velocities of origin and target orbit within their system
    origin_v_esc = vesc_r(origin_mu, origin_radius)
    target_v_esc = vesc_r(target_mu, target_radius)

    # Velocity in the peri-center of the parabolic orbits within their system
    v2p = v_vesc(origin_v_esc, hyperbolic_v_inf_origin)
    v4p = v_vesc(target_v_esc, hyperbolic_v_inf_target)

    # Calculate the velocities of the origin and target orbits within their system
    v1 = orbit_velocity(origin_mu, origin_radius, origin_a)
    v5 = orbit_velocity(target_mu, target_radius, target_a)

    # Calculate the required velocity differences in the peri-center of origin and target orbit within their systems
    dv12p = delta(v1, v2p)
    dv4p5 = delta(v4p, v5)

    return dv12p, dv4p5
