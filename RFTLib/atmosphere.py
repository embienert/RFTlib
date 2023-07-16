import numpy as np
import math

from .util import EARTH_RADIUS, EARTH_MU


M = 0.028_964_4  # Molar mass of dry air in kg/mol
R = 8.31432  # Universal gas coefficient in J/(mol*K)
kappa = 1.4  # ideal two-atom isoentrope exponent


# Atmosphere default altitude and temperature values
TROPOSPHERE_ALTITUDE = 0  # m
TROPOSPHERE_STANDARD_T = 288.15  # K
TROPOSPHERE_SINK_RATE = -0.0065  # K/m

TROPOPAUSE_ALTITUDE = 11_000  # m
TROPOPAUSE_STANDARD_T = 216.65  # K
TROPOPAUSE_SINK_RATE = 0.0  # K/m

STRATOSPHERE_ALTITUDE = 20_000  # m
STRATOSPHERE_STANDARD_T = 216.65  # K
STRATOSPHERE_SINK_RATE = 0.001  # K/m

STRATOSPHERE2_ALTITUDE = 32_000  # m
STRATOSPHERE2_STANDARD_T = 228.65  # K
STRATOSPHERE2_SINK_RATE = 0.0028  # K/m

STRATOPAUSE_ALTITUDE = 47_000  # m
STRATOPAUSE_STANDARD_T = 270.65  # K
STRATOPAUSE_SINK_RATE = 0.0  # K/m

MESOSPHERE_ALTITUDE = 51_000  # m
MESOSPHERE_STANDARD_T = 270.65  # K
MESOSPHERE_SINK_RATE = -0.0028  # K/m

MESOSPHERE2_ALTITUDE = 71_000  # m
MESOSPHERE2_STANDARD_T = 214.65  # K
MESOSPHERE2_SINK_RATE = -0.002  # K/m


def local_gravity(altitude: float):
    """
    Calculate the local gravity in earth's influence

    :param altitude: [m] Altitude above earth's sea level
    :return: g [m/s²] Gravity acceleration
    """

    # Calculation
    local_radius = EARTH_RADIUS + altitude

    return EARTH_MU / local_radius ** 2


def local_atmospheric_temperature(altitude: float):
    """
    Calculate the local atmospheric temperature

    :param altitude: [m] Altitude above sea level
    :return: T [K] Local atmospheric temperature
    """

    LIMIT = 80_000  # m


    # Calculations
    if altitude < TROPOPAUSE_ALTITUDE:
        # Troposphere
        return TROPOSPHERE_STANDARD_T + (altitude - TROPOSPHERE_ALTITUDE) * TROPOSPHERE_SINK_RATE
    elif altitude < STRATOSPHERE_ALTITUDE:
        # Tropopause
        return TROPOPAUSE_STANDARD_T + (altitude - TROPOPAUSE_ALTITUDE) * TROPOPAUSE_SINK_RATE
    elif altitude < STRATOSPHERE2_ALTITUDE:
        # Stratosphere
        return STRATOSPHERE_STANDARD_T + (altitude - STRATOSPHERE_ALTITUDE) * STRATOSPHERE_SINK_RATE
    elif altitude < STRATOPAUSE_ALTITUDE:
        # Stratosphere 2
        return STRATOSPHERE2_STANDARD_T + (altitude - STRATOSPHERE2_ALTITUDE) * STRATOSPHERE2_SINK_RATE
    elif altitude < MESOSPHERE_ALTITUDE:
        # Stratopause
        return STRATOPAUSE_STANDARD_T + (altitude - STRATOPAUSE_ALTITUDE) * STRATOPAUSE_SINK_RATE
    elif altitude < MESOSPHERE2_ALTITUDE:
        # Mesosphere
        return MESOSPHERE_STANDARD_T + (altitude - MESOSPHERE_ALTITUDE) * MESOSPHERE_SINK_RATE
    elif altitude <= LIMIT:
        # Mesosphere 2
        return MESOSPHERE2_STANDARD_T + (altitude - MESOSPHERE2_ALTITUDE) * MESOSPHERE2_SINK_RATE
    else:
        # Not valid above LIMIT altitude
        raise ValueError(f"Altitude must be at most {LIMIT}m")


def density_no_dt(rho_b: float, h_b: float, altitude: float):
    """
    Calculate the local atmospheric density, assuming a temperature sink rate of 0 K/m

    :param rho_b: [kg/m³] Reference barometric air density
    :param h_b: [m] Reference altitude above sea level
    :param altitude: [m] Altitude above sea level
    :return: rho [kg/m³] Barometric air density
    """

    T_b = local_atmospheric_temperature(h_b)
    g0 = local_gravity(altitude)

    return rho_b * np.exp(-g0 * M * (altitude - h_b) / (R * T_b))


def density_dt(rho_b: float, dt: float, h_b: float, altitude: float):
    """
    Calculate the local atmospheric density, assuming a non-zero temperature sink rate dt

    :param rho_b: [kg/m³] Reference barometric air density
    :param dt: [K/m] Temperature sink date
    :param h_b: [m] Reference altitude above sea level
    :param altitude: [m] Altitude above sea level
    :return: rho [kg/m³] Barometric air density
    """

    T_b = local_atmospheric_temperature(h_b)
    T = local_atmospheric_temperature(altitude)
    g0 = local_gravity(altitude)


    return rho_b * math.pow(T_b / T, 1 + g0 * M / (R * dt))


def density(rho_b: float, dt: float, h_b: float, altitude: float):
    """
    Calculate the local atmospheric density

    :param rho_b: [kg/m³] Reference barometric air density
    :param dt: [K/m] Temperature sink date
    :param h_b: [m] Reference altitude above sea level
    :param altitude: [m] Altitude above sea level
    :return: rho [kg/m³] Barometric air density
    """

    if dt == 0:
        return density_no_dt(rho_b, h_b, altitude)
    else:
        return density_dt(rho_b, dt, h_b, altitude)


def local_atmospheric_density(altitude: float):
    """
    Calculate the local atmospheric density

    :param altitude: [m] Altitude above sea level
    :return: rho [kg/m³] Barometric air density
    """

    # Constants
    TROPOSPHERE_ALTITUDE = 0  # m
    TROPOSPHERE_SINK_RATE = -0.0065  # K/m
    TROPOSPHERE_STANDARD_DENSITY = 1.22500  # kg/m³

    TROPOPAUSE_ALTITUDE = 11_000  # m
    TROPOPAUSE_SINK_RATE = 0.0  # K/m
    TROPOPAUSE_STANDARD_DENSITY = 0.363910  # kg/m³

    STRATOSPHERE_ALTITUDE = 20_000  # m
    STRATOSPHERE_SINK_RATE = 0.001  # K/m
    STRATOSPHERE_STANDARD_DENSITY = 0.088030  # kg/m³

    STRATOSPHERE2_ALTITUDE = 32_000  # m
    STRATOSPHERE2_SINK_RATE = 0.0028  # K/m
    STRATOSPHERE2_STANDARD_DENSITY = 0.013220  # kg/m³

    STRATOPAUSE_ALTITUDE = 47_000  # m
    STRATOPAUSE_SINK_RATE = 0.0  # K/m
    STRATOPAUSE_STANDARD_DENSITY = 0.001430  # kg/m³

    MESOSPHERE_ALTITUDE = 51_000  # m
    MESOSPHERE_SINK_RATE = -0.0028  # K/m
    MESOSPHERE_STANDARD_DENSITY = 0.000860  # kg/m³

    MESOSPHERE2_ALTITUDE = 71_000  # m
    MESOSPHERE2_SINK_RATE = -0.002  # K/m
    MESOSPHERE2_STANDARD_DENSITY = 0.000064  # kg/m³

    LIMIT = 80_000  # m


    # Calculations
    if altitude < TROPOPAUSE_ALTITUDE:
        # Troposphere
        return density(TROPOSPHERE_STANDARD_DENSITY, TROPOSPHERE_SINK_RATE, TROPOSPHERE_ALTITUDE, altitude)
    elif altitude < STRATOSPHERE_ALTITUDE:
        # Tropopause
        return density(TROPOPAUSE_STANDARD_DENSITY, TROPOPAUSE_SINK_RATE, TROPOPAUSE_ALTITUDE, altitude)
    elif altitude < STRATOSPHERE2_ALTITUDE:
        # Stratosphere
        return density(STRATOSPHERE_STANDARD_DENSITY, STRATOSPHERE_SINK_RATE, STRATOSPHERE_ALTITUDE, altitude)
    elif altitude < STRATOPAUSE_ALTITUDE:
        # Stratosphere 2
        return density(STRATOSPHERE2_STANDARD_DENSITY, STRATOSPHERE2_SINK_RATE, STRATOSPHERE2_ALTITUDE, altitude)
    elif altitude < MESOSPHERE_ALTITUDE:
        # Stratopause
        return density(STRATOPAUSE_STANDARD_DENSITY, STRATOPAUSE_SINK_RATE, STRATOPAUSE_ALTITUDE, altitude)
    elif altitude < MESOSPHERE2_ALTITUDE:
        # Mesosphere
        return density(MESOSPHERE_STANDARD_DENSITY, MESOSPHERE_SINK_RATE, MESOSPHERE_ALTITUDE, altitude)
    elif altitude <= LIMIT:
        # Mesosphere 2
        return density(MESOSPHERE2_STANDARD_DENSITY, MESOSPHERE2_SINK_RATE, MESOSPHERE2_ALTITUDE, altitude)
    else:
        # Not valid above LIMIT altitude
        raise ValueError(f"Altitude must be at most {LIMIT}m")


def local_atmospheric_pressure(altitude: float):
    """
    Calculate the local atmospheric pressure

    :param altitude: [m] Altitude above sea level
    :return: p [Pa] Barometric pressure
    """

    rho = local_atmospheric_density(altitude)
    T = local_atmospheric_temperature(altitude)

    return rho * R * T / M


def local_sonic_speed(altitude: float):
    """
    Calculate the local sonic velocity

    :param altitude: [m] Altitude above sea level
    :return: c [m/s] Local sonic velocity
    """

    pressure = local_atmospheric_pressure(altitude)
    rho = local_atmospheric_density(altitude)

    return np.sqrt(kappa * pressure / rho)
