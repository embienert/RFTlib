import numpy as np


def linear_eccentricity_a_b(a: float, b: float) -> float:
    """
    Calculate the ellipse's linear eccentricity

    :param a: [km] Semi-major axis of the ellipse
    :param b: [km] Semi-minor axis of the ellipse
    :return: c [km] Linear eccentricity of the ellipse
    """

    return np.sqrt(a**2 - b**2)


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
