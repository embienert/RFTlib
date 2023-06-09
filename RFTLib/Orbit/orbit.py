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


def orbit_period_mu(mu: float, a: float) -> float:
    """
    Calculate the orbit period of an object around a reference object, assuming the object's orbit is circular (r=a)

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param a: [km] Semi-major axis of the orbit around the reference object
    :return: T [s] Orbit period of the object
    """

    return 2 * np.pi * np.sqrt(a**3 / mu)




def orbit_period_circular_v(r: float, v: float) -> float:
    """
    Calculate the orbit period of an object around a reference object, assuming the object's orbit is circular (r=a)

    :param r: [km] Radius of the orbit around the reference object
    :param v: [km/s] Velocity of the orbiting object
    :return: T [s] Orbit period of the object
    """

    return 2 * np.pi * r / v


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


def true_anomaly(rx: float, ry: float) -> float:
    """
    Calculate the true anomaly of an orbiting object

    :param rx: [km] Current offset of the orbiting object along the Euclidean x-Axis on the orbit plane
    :param ry: [km] Current offset of the orbiting object along the Euclidean y-Axis on the orbit plane
    :return: nu [rad] True anomaly of the orbiting object
    """

    return np.arctan(ry / rx)


def aop(rx_p: float, ry_p: float) -> float:
    """
    Calculate the argument of perigee (aop) of an object's orbit

    :param rx_p: [km] Current offset of the orbiting object along the Euclidean x-Axis on the orbit plane
    :param ry_p: [km] Current offset of the orbiting object along the Euclidean y-Axis on the orbit plane
    :return: omega [rad] Argument of perigee (aop) of the orbit
    """

    return np.arctan(ry_p / rx_p)


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


def vesc_r(mu: float, r: float) -> float:
    """
    Calculate velocity required to exit the gravitational influence of a reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param r: [km] Radius of the orbit around the reference object
    :return: v_esc [km/s] Escape velocity required to exit SOI of the reference object
    """

    return np.sqrt(2 * mu / r)


def v_vesc(v_esc: float, v_inf: float) -> float:
    """
    Calculate required velocity to achieve the desired hyperbolic excess speed

    :param v_esc: [km/s] Escape velocity required to exit SOI of the reference object
    :param v_inf: [km/s] Hyperbolic excess velocity
    :return: v [km/s] Required velocity to achieve desired hyperbolic excess velocity
    """

    return np.sqrt(v_esc**2 + v_inf**2)


def r_SOI_mu(a: float, mu: float, MU: float) -> float:
    """
    Calculate the radius of the sphere of influence (SOI) of an object orbiting around a reference object

    :param a: [km] Semi-major axis of the object's orbit around the reference object
    :param mu: [km³/s²] Standard gravitational parameter of the orbiting object
    :param MU: [km³/s²] Standard gravitational parameter of the orbited reference object
    :return: r [km] Radius of the sphere of influence of the orbiting object
    """

    return a * np.pow(mu/MU, 2/5)


def r_SOI_m(a: float, m: float, M: float) -> float:
    """
    Calculate the radius of the sphere of influence (SOI) of an object orbiting around a reference object

    :param a: [km] Semi-major axis of the object's orbit around the reference object
    :param m: [kg] Mass of the orbiting object
    :param M: [kg] Mass of the orbited reference object
    :return: r [km] Radius of the sphere of influence of the orbiting object
    """

    return a * np.pow(m/M, 2/5)


def energy_v_m_r(MU: float, v: float, m: float, r: float) -> float:
    """
    Calculate the energy of an object in the gravitational field of a reference object

    :param MU: [km³/s²] Standard gravitational parameter of the reference object
    :param v: [km/s] Velocity of the object
    :param m: [kg] Mass of the object
    :param r: [km] Distance of the object from the reference object
    :return: E [kg*km²/s² bzw. MJ] Orbital energy of the object
    """

    return (m * v**2 / 2) - (m * MU / r)


def sam(mu: float, p: float) -> float:
    """
    Calculate the specific angular momentum (sam) of an object orbiting a reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param p: [km] Semilatus rectum of the object's orbit
    :return: h [km²/s] Specific angular momentum of the object
    """

    return np.sqrt(mu * p)


def p_sam(mu: float, h: float) -> float:
    """
    Calculate the semilatus rectum of an object's orbit around a reference object

    :param mu: [km³/s²] Standard gravitational parameter of the reference object
    :param h: [km²/s] Specific angular momentum (sam) of the orbiting object
    :return: p [km] Semilatus rectum of the object's orbit
    """

    return h**2 / mu


def mu_sam_p(h: float, p: float) -> float:
    """
    Calculate a reference object's standard gravitational parameter by an orbiting object's attributes

    :param h: [km²/s] Specific angular momentum (sam) of the orbiting object
    :param p: [km] Semilatus rectum of the object's orbit
    :return: mu [km³/s²] Standard gravitational parameter of the reference object
    """

    return h**2 / p


def soe(MU: float, v: float, r: float):
    """
    Calculate the specific orbital energy of an object orbiting a reference object

    :param MU: [km³/s²] Standard gravitational parameter of the reference object
    :param v: [km/s] Velocity of the object
    :param r: [km] Distance of the object from the reference object
    :return: e [km²/s² bzw. MJ/kg] Specific orbital energy of the object
    """

    return (v**2 / 2) - (MU / r)


def soe_E_m(E: float, m: float) -> float:
    """
    Calculate the specific orbital energy of an object orbiting a reference object

    :param E: [kg*km²/s² bzw. MJ] Orbital energy of the orbiting object
    :param m: [kg] Mass of the orbiting object
    :return: e [km²/s² bzw. MJ/kg] Specific orbital energy of the object
    """

    return E / m


def E_soe_m(e: float, m: float) -> float:
    """
    Calculate the orbital energy of an object orbiting around a reference object

    :param e: [km²/s² bzw MJ/kg] Specific orbital energy of the object
    :param m: [kg] Mass of the orbiting object
    :return: E [kg*km²/s² bzw. MJ] Orbital energy of the orbiting object
    """

    return e * m


def m_soe_E(e: float, E: float) -> float:
    """
    Calculate the mass of an object orbiting a reference object

    :param e: [km²/s² bzw MJ/kg] Specific orbital energy of the object
    :param E: [kg*km²/s² bzw. MJ] Orbital energy of the orbiting object
    :return: [kg] Mass of the orbiting object
    """

    return E / e


def sam_SAM_m(H: float, m: float) -> float:
    """
    Calculate the specific angular momentum of an object orbiting a reference object

    :param H: [kg*km²/s bzw. kNms] constant angular momentum of the reference object
    :param m: [kg] Mass of the object
    :return: h [km²/s] Specific angular momentum of the object
    """

    return H / m


def SAM_sam_m(h: float, m: float) -> float:
    """
    Calculate the constant angular momentum of a reference object orbited by another object

    :param h: [km²/s] Specific angular momentum of the object
    :param m: [kg] Mass of the object
    :return: H [kg*km²/s bzw. kNms] constant angular momentum of the reference object
    """

    return h * m


def m_SAM_sam(H: float, h: float) -> float:
    """
    Calculate the mass of an object orbiting a reference object

    :param H: [kg*km²/s bzw. kNms] constant angular momentum of the reference object
    :param h: [km²/s] Specific angular momentum of the object
    :return: [kg] Mass of the object
    """

    return H / h


def sam_rp_vp(radius: float, velocity: float) -> float:
    """
    Calculate the specific angular momentum of an object orbiting a reference object

    :param radius: [km] Distance of the object in its peri- or apo-center
    :param velocity: [km/s] Velocity of the object in its peri- or apo-center
    :return: [km²/s] Specific angular momentum of the object
    """

    return radius * velocity
