import numpy as np


def linear_eccentricity_a_b(a: float, b: float) -> float:
    """
    Calculate the ellipse's linear eccentricity

    :param a: [km] Semi-major axis of the ellipse
    :param b: [km] Semi-minor axis of the ellipse
    :return: c [km] Linear eccentricity of the ellipse
    """

    return np.sqrt(a**2 - b**2)


def semimajor_axis_b_c(b: float, c: float) -> float:
    """
    Calculate the ellipse's semi-major axis

    :param b: [km] Semi-minor axis of the ellipse
    :param c: [km] Linear eccentricity of the ellipse
    :return: a [km] Semi-major axis of the ellipse
    """

    return np.sqrt(c**2 + b**2)


def semiminor_axis_a_c(a: float, c: float) -> float:
    """
    Calculate the ellipse's semi-minor axis

    :param a: [km] Semi-major axis of the ellipse
    :param c: [km] Linear eccentricity of the ellipse
    :return: b [km] Semi-minor axis of the ellipse
    """

    return np.sqrt(a**2 - c**2)


def numeric_eccentricity_a_c(a: float, c: float) -> float:
    """
    Calculate the ellipse's numeric eccentricity

    NOTE:
        e=0 -> Circle
        e=1 -> Parabola

    :param a: [km] Semi-major axis of the ellipse
    :param c: [km] Linear eccentricity of the ellipse
    :return: e [] Numeric eccentricity of the ellipse (between 0 and 1)
    """

    return c / a


def numeric_eccentricity_rp_ra(radius_peri: float, radius_apo: float) -> float:
    """
    Calculate the numeric eccentricity

    :param radius_peri: [km] Radius of the ellipse's peri-center
    :param radius_apo: [km] Radius of the ellipse's apo-center
    :return: e [km] Linear eccentricity of the ellipse
    """

    a = (radius_peri + radius_apo) / 2
    return (radius_apo - a) / a


def linear_eccentricity_rp_ra(radius_peri: float, radius_apo: float) -> float:
    """
    Calculate the numeric eccentricity

    :param radius_peri: [km] Radius of the ellipse's peri-center
    :param radius_apo: [km] Radius of the ellipse's apo-center
    :return: c [km] Linear eccentricity of the ellipse
    """

    return radius_apo - (radius_peri + radius_apo) / 2


def semilatus_rectum_a_b(a: float, b: float) -> float:
    """
    Calculate the ellipse's semilatus rectum ("Halbparameter")

    :param a: [km] Semi-major axis of the ellipse
    :param b: [km] Semi-minor axis of the ellipse
    :return: p [km] Semilatus rectum of the ellipse
    """

    return b**2 / a


def semilatus_rectum_a_e(a: float, e: float) -> float:
    """
    Calculate the ellipse's semilatus rectum ("Halbparameter")

    :param a: [km] Semi-major axis of the ellipse
    :param e: [] Numeric eccentricity of the ellipse
    :return: p [km] Semilatus rectum of the ellipse
    """

    return a * (1 - e**2)


def radius_peri_a_e(a: float, e: float) -> float:
    """
    Calculate the radius of the ellipse in its peri-center

    :param a: [km] Semi-major axis of the ellipse
    :param e: [] Numeric eccentricity of the ellipse
    :return: r_p [km] Radius of the ellipse in its peri-center
    """

    return a * (1 - e)


def radius_apo_a_e(a: float, e: float) -> float:
    """
    Calculate the radius of the ellipse in its apo-center

    :param a: [km] Semi-major axis of the ellipse
    :param e: [] Numeric eccentricity of the ellipse
    :return: r_a [km] Radius of the ellipse in its apo-center
    """

    return a * (1 + e)


def orbit_period_kepler3(a1: float, a2: float, T2: float) -> float:
    """
    Calculate the orbit time using Kepler's third law

    :param a1: [km] Semi-major axis of the orbit with unknown orbit period
    :param a2: [km] Semi-major axis of the orbit with known orbit period
    :param T2: [s] Orbit period of the orbit corresponding to a2
    :return: T1 [s] Orbit period of the orbit corresponding to a1
    """

    return np.sqrt((a1 / a2)**3 * T2**2)


def a_kepler3(a2: float, T1: float, T2: float) -> float:
    """
    Calculate the semi-major axis using Kepler's third law

    :param a2: [km] Semi-major axis of the orbit corresponding to T2
    :param T1: [km] Orbit period of the orbit with unknown semi-major axis
    :param T2: [km] Orbit period of the orbit with known semi-major axis
    :return: a1 [km] Semi-major axis of the orbit corresponding to T1
    """

    return float(np.cbrt((T1 / T2)**2 * a2**3))
