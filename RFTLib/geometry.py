import numpy as np


def sphere_surface(radius: float) -> float:
    """
    Calculates a sphere's surface

    :param radius: [km] Sphere radius
    :return: A [km²] Sphere surface area
    """

    return 4 * np.pi * radius**2


def sphere_volume(radius: float) -> float:
    """
    Calculates a sphere's volume

    :param radius: [km] Sphere radius
    :return: V [km³] Sphere volume
    """

    return 4/3 * np.pi * radius**3


def sphere_segment_surface_theta(radius: float, theta: float) -> float:
    """
    Calculates the mantle surface of a sphere cap

    :param radius: [km] Sphere radius
    :param theta: [rad] Polar angle
    :return: A [km²] Mantle surface area of the spherical cap
    """

    return np.pi * radius**2 * (1 - np.cos(theta))


def sphere_segment_surface_base(height: float, base_radius: float) -> float:
    """
    Calculates the mantle surface of a sphere cap

    :param height: [km] Height of the spherical cap (not distance from center of the sphere!)
    :param base_radius: [km] Radius of the base plane of the spherical cap
    :return: A [km²] Mantle surface area of the spherical cap
    """

    return np.pi * (height**2 + base_radius**2)


def sphere_segment_surface_height(radius: float, height: float) -> float:
    """
    Calculates the mantle surface of a sphere cap

    :param radius: [km] Sphere radius
    :param height: [km] Distance from the center of the sphere to the sphere cap base
    :return: A [km²] Mantle surface area of the spherical cap
    """

    return 2 * np.pi * radius * (radius - height)
