from typing import List

import numpy as np


def dv_veff_m(v_eff: float, m_start: float, m_end: float) -> float:
    """
    Calculate the rocket's velocity difference

    :param v_eff: [km/s] Effective exit velocity of the propulsion system
    :param m_start: [kg] Initial mass of the rocket
    :param m_end: [kg] Final mass of the rocket
    :return: dv [km/s] Velocity difference of the rocket after engine burn
    """

    return v_eff * np.log(m_start / m_end)


def veff_dv_m(dv: float, m_start: float, m_end: float) -> float:
    """
    Calculate the propulsion system's effective exit velocity

    :param dv: [km/s] Velocity difference of the rocket after engine burn
    :param m_start: [kg] Initial mass of the rocket
    :param m_end: [kg] Final mass of the rocket
    :return: v_eff [km/s] Effective exit velocity of the propulsion system
    """

    return dv / np.log(m_start / m_end)


def dv_veff_mpf(v_eff: float, mpf: float) -> float:
    """
    Calculate the rocket's velocity difference

    :param v_eff: [km/s] Effective exit velocity of the propulsion system
    :param mpf: [] Propellant mass fraction of the burn process
    :return: dv [km/s] Velocity difference of the rocket after engine burn
    """

    return v_eff * np.log(1 / (1 - mpf))


def mpf_dv_veff(dv: float, veff: float) -> float:
    """
    Calculate propellant mass fraction

    :param dv: [km/s] Velocity difference of the rocket after engine burn
    :param veff: [km/s] Effective exit velocity of the propulsion system
    :return: mpf [] Propellant mass fraction of the burn process
    """

    return 1 - np.power(np.e, -dv / veff)


def veff_dv_mpf(dv: float, mpf: float) -> float:
    """
    Calculate the propulsion system's effective exit velocity

    :param dv: [km/s] Velocity difference of the rocket after engine burn
    :param mpf: [] Propellant mass fraction of the burn process
    :return: v_eff [km/s] Effective exit velocity of the propulsion system
    """

    return dv / np.log(1 / (1 - mpf))


def mpf_fuel_mass_m0(fuel_mass: float, m0: float) -> float:
    """
    Calculate propellant mass fraction

    :param fuel_mass: [kg] Propellant mass that is burned
    :param m0: [kg] Initial mass of the rocket (including propellant mass)
    :return: mpf [] Propellant mass fraction of the burn process
    """

    return fuel_mass / m0


def fuel_mass_mpf_m0(mpf: float, m0: float) -> float:
    """
    Calculate the propellant mass that is burned

    :param mpf: [] Propellant mass fraction of the burn process
    :param m0: [kg] Initial mass of the rocket (including propellant mass)
    :return: m_f [kg] Propellant mass that is burned
    """

    return mpf * m0


def m0_mpf_fuel_mass(mpf: float, fuel_mass: float) -> float:
    """
    Calculate the initial mass of the rocket

    :param mpf: [] Propellant mass fraction of the burn process
    :param fuel_mass: [kg] Propellant mass that is burned
    :return: m0 [kg] Initial mass of the rocket (including propellant mass)
    """

    return fuel_mass / mpf


def veff_isp(g0: float, isp: float) -> float:
    """
    Calculate the propulsion system's effective exit velocity

    :param g0: [km/s²] Standard acceleration of gravity
    :param isp: [s] Specific impulse of the propulsion system
    :return: veff [km/s] Effective exit velocity of the propulsion system
    """

    return g0 * isp


def isp_veff(g0: float, v_eff: float) -> float:
    """
    Calculate the specific impulse of the propulsion system

    :param g0: [km/s²] Standard acceleration of gravity
    :param v_eff: [km/s] Effective exit velocity of the propulsion system
    :return: isp [s] Specific impulse of the propulsion system
    """

    return v_eff / g0


def veff_thrust_flow_rate(thrust: float, flow_rate: float) -> float:
    """
    Calculate the propulsion system's effective exit velocity

    :param thrust: [km*kg/s² bzw. kN] Thrust of the propulsion system
    :param flow_rate: [kg/s] Propellant flow rate of the propulsion system
    :return: veff [km/s] Effective exit velocity of the propulsion system
    """

    return thrust / flow_rate


def isp_thrust_flow_rate(g0: float, thrust: float, flow_rate: float) -> float:
    """
    Calculate the specific impulse of the propulsion system

    :param g0: [km/s²] Standard acceleration of gravity
    :param thrust: [km*kg/s² bzw. kN] Thrust of the propulsion system
    :param flow_rate: [kg/s] Propellant flow rate of the propulsion system
    :return: isp [s] Specific impulse of the propulsion system
    """

    return thrust / (g0 * flow_rate)


def flow_rate_isp_thrust(g0: float, isp: float, thrust: float) -> float:
    """
    Calculate the propellant flow rate of the propulsion system

    :param g0: [km/s²] Standard acceleration of gravity
    :param isp: [s] Specific impulse of the propulsion system
    :param thrust: [km*kg/s² bzw. kN] Thrust of the propulsion system
    :return: ṁ [kg/s] Propellant flow rate of the propulsion system
    """

    return thrust / (g0 * isp)


def thrust_isp_flow_rate(v_eff: float, flow_rate: float) -> float:
    """
    Calculate the thrust of the propulsion system

    :param v_eff: [km/s] Effective exit velocity of the propulsion system
    :param flow_rate: [kg/s] Propellant flow rate of the propulsion system
    :return: T [km*kg/s² bzw. kN] Thrust of the propulsion system
    """

    return v_eff * flow_rate


# TODO: Calculating total thrust/effective velocity from multiple parallel propulsion systems
def veff_total(v_effs: List[float], flow_rates: List[float]) -> float:
    """
    Calculate the total effective exit velocity of a parallel propulsion system

    :param v_effs: *[km/s] List containing the effective exit velocities of the single propulsion systems
    :param flow_rates: *[kg/s] List containing the propellant flow rates of the single propulsion systems
    :return: v_eff [km/s] Total effective exit velocity of all parallel propulsion systems
    """

    assert len(v_effs) == len(flow_rates), "List sizes must be equal"

    v_effs_array = np.array(v_effs)
    flow_rates_array = np.array(flow_rates)

    return np.sum(v_effs_array * flow_rates_array) / np.sum(flow_rates_array)


def veff_veff_total(v_eff_total: float, v_effs: List[float], flow_rates: List[float]) -> float:
    """
    Calculate the effective exit velocity of one propulsion system in a set of parallel propulsion systems

    :param v_eff_total: [km/s] Total effective exit velocity of the parallel propulsion system
    :param v_effs: *[km/s] List containing the effective exit velocities of the single propulsion systems
        (except for the wanted value)
    :param flow_rates: *[kg/s] List containing the propellant flow rates of the single propulsion systems,
        with the first entry corresponding to the wanted effective exit velocity
    :return: v_eff [km/s] Effective exit velocity of the single propulsion system
    """

    assert len(flow_rates) - len(v_effs) == 1, \
        "There must be exactly one more entry in the flow_rates list than the v_effs list"

    v_effs_array = np.array(v_effs)
    flow_rates_array = np.array(flow_rates)

    return (v_eff_total * np.sum(flow_rates_array) - np.sum(v_effs_array * flow_rates_array[1:])) / flow_rates_array[0]


def flow_rate_veff_total(v_eff_total: float, v_effs: List[float], flow_rates: List[float]) -> float:
    """
    Calculate the propellant flow rate of one propulsion system in a set of parallel propulsion systems

    :param v_eff_total: [km/s] Total effective exit velocity of the parallel propulsion system
    :param v_effs: *[km/s] List containing the effective exit velocities of the single propulsion systems,
        with the first entry corresponding to the wanted flow rate
    :param flow_rates: *[kg/s] List containing the propellant flow rates of the single propulsion systems
        (except for the wanted value)
    :return: ṁ [kg/s] Propellant flow rate of the single propulsion system
    """

    assert len(v_effs) - len(flow_rates) == 1, \
        "There must be exactly one more entry in the v_effs list than the flow_rates list"

    v_effs_array = np.array(v_effs)
    flow_rates_array = np.array(flow_rates)

    return (v_eff_total * np.sum(flow_rates_array) - np.sum(flow_rates_array * v_effs_array[1:])) / \
        (v_effs[0] - v_eff_total)


def structure_mass_coefficient(structure_mass: float, fuel_mass: float) -> float:
    """
    Calculate the structure mass coefficient

    :param structure_mass: [kg] Structure and motor mass of the rocket
    :param fuel_mass: [kg] Propellant mass of the rocket
    :return: sigma [] Structure mass coefficient of the rocket
    """

    return structure_mass / (structure_mass + fuel_mass)


def structure_mass_sigma(sigma: float, fuel_mass: float) -> float:
    """
    Calculate the structure mass

    :param sigma: [] Structure mass coefficient of the rocket
    :param fuel_mass: [kg] Propellant mass of the rocket
    :return: m_S [kg] Structure and motor mass of the rocket
    """

    return sigma * fuel_mass / (1 - sigma)


def fuel_mass_sigma(sigma: float, structure_mass: float) -> float:
    """
    Calculate the fuel mass

    :param sigma: [] Structure mass coefficient of the rocket
    :param structure_mass: [kg] Structure and motor mass of the rocket
    :return: m_T [kg] Propellant mass of the rocket
    """

    return (1 - sigma) * structure_mass / sigma


def itot_thrust_burn_duration(thrust: float, burn_duration: float) -> float:
    """
    Calculate the total impulse

    :param thrust: [kg*km/s² bzw. kN] Thrust of the propulsion system
    :param burn_duration: [s] Duration of the burn process
    :return: I_tot [kg*km/s bzw. kNs] Total impulse of the propulsion system
    """

    return thrust * burn_duration


def thrust_itot_burn_duration(i_tot: float, burn_duration: float) -> float:
    """
    Calculate the thrust of the propulsion system

    :param i_tot: [kg*km/s bzw. kNs] Total impulse of the propulsion system
    :param burn_duration: [s] Duration of the burn process
    :return: T [kg*km/s² bzw. kN] Thrust of the propulsion system
    """

    return i_tot / burn_duration


def burn_duration_itot_thrust(i_tot: float, thrust: float) -> float:
    """
    Calculate the burn time of the propulsion system

    :param i_tot: [kg*km/s bzw. kNs] Total impulse of the propulsion system
    :param thrust: [kg*km/s² bzw. kN] Thrust of the propulsion system
    :return: t [s] Duration of the burn process
    """

    return i_tot / thrust


def itot_veff_m(v_eff: float, m: float) -> float:
    """
    Calculate the total impulse

    :param v_eff: [km/s] Effective exit velocity of the propulsion system
    :param m: [kg] Mass of the fuel used during the burn process
    :return: I_tot [kg*km/s bzw. kNs] Total impulse of the propulsion system
    """

    return v_eff * m


def veff_itot_m(i_tot: float, m: float) -> float:
    """
    Calculate the effective exit velocity of the propulsion system

    :param i_tot: [kg*km/s bzw. kNs] Total impulse of the propulsion system
    :param m: [kg] Mass of the fuel used during the burn process
    :return: v_eff [km/s] Effective exit velocity of the propulsion system
    """

    return i_tot / m


def m_itot_veff(i_tot: float, v_eff: float) -> float:
    """
    Calculate the mass of the fuel used during the burn process

    :param i_tot: [kg*km/s bzw. kNs] Total impulse of the propulsion system
    :param v_eff: [km/s] Effective exit velocity of the propulsion system
    :return: m [kg] Mass of the fuel used during the burn process
    """

    return i_tot / v_eff


def veff_saint_venant_wantzel(kappa: float, R: float, M: float, T_c: float, p_e: float, p_c: float) -> float:
    """
    Calculate the effective exit velocity of a propulsion system according to the equation of Saint-Venant and Wantzel

    :param kappa: [] Adiabatic coefficient of the propellant
    :param R: [kg*km²/(s²*K*mol)] Universal gas constant
    :param M: [kg/mol] Molar mass of the propellant
    :param T_c: [K] Temperature in the burn chamber
    :param p_e: [kg/(km*s²)] Pressure at the exhaust
    :param p_c: [kg/(km*s²)] Pressure in the burn chamber
    :return: v_eff [km/s] Effective exit velocity of the propulsion system
    """

    return np.sqrt(2 * kappa * R * T_c / ((kappa - 1) * M) * (1 - np.pow(p_e / p_c, (kappa - 1) / kappa)))


def v_launch_latitude(rotation_velocity: float, latitude: float) -> float:
    """
    Calculate the velocity gain when launching east.

    :param rotation_velocity: [km/s] Rotation velocity of the body the rocket is launching from
    :param latitude: [rad] Latitude of the position of the rocket launch
    :return: v [km/s] Velocity gain from launching east
    """

    return rotation_velocity * np.cos(latitude)
