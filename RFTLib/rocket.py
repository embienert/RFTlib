from typing import List

import numpy as np


def dv_veff_m(v_eff: float, m_start: float, m_end: float) -> float:
    """
    Calculate the rocket's velocity difference

    :param v_eff: [m/s] Effective exit velocity of the propulsion system
    :param m_start: [kg] Initial mass of the rocket
    :param m_end: [kg] Final mass of the rocket
    :return: dv [m/s] Velocity difference of the rocket after engine burn
    """

    return v_eff * np.log(m_start / m_end)


def veff_dv_m(dv: float, m_start: float, m_end: float) -> float:
    """
    Calculate the propulsion system's effective exit velocity

    :param dv: [m/s] Velocity difference of the rocket after engine burn
    :param m_start: [kg] Initial mass of the rocket
    :param m_end: [kg] Final mass of the rocket
    :return: v_eff [m/s] Effective exit velocity of the propulsion system
    """

    return dv / np.log(m_start / m_end)


def dv_veff_mpf(v_eff: float, mpf: float) -> float:
    """
    Calculate the rocket's velocity difference

    :param v_eff: [m/s] Effective exit velocity of the propulsion system
    :param mpf: [] Propellant mass fraction of the burn process
    :return: dv [m/s] Velocity difference of the rocket after engine burn
    """

    return v_eff * np.log(1 / (1 - mpf))


def mpf_dv_veff(dv: float, veff: float) -> float:
    """
    Calculate propellant mass fraction

    :param dv: [m/s] Velocity difference of the rocket after engine burn
    :param veff: [m/s] Effective exit velocity of the propulsion system
    :return: mpf [] Propellant mass fraction of the burn process
    """

    return 1 - np.power(np.e, -dv/veff)


def veff_dv_mpf(dv: float, mpf: float) -> float:
    """
    Calculate the propulsion system's effective exit velocity

    :param dv: [m/s] Velocity difference of the rocket after engine burn
    :param mpf: [] Propellant mass fraction of the burn process
    :return: v_eff [m/s] Effective exit velocity of the propulsion system
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

    :param g0: [m/s²] Standard acceleration of gravity
    :param isp: [s] Specific impulse of the propulsion system
    :return: veff [m/s] Effective exit velocity of the propulsion system
    """

    return g0 * isp


def isp_veff(g0: float, v_eff: float) -> float:
    """
    Calculate the specific impulse of the propulsion system

    :param g0: [m/s²] Standard acceleration of gravity
    :param v_eff: [m/s] Effective exit velocity of the propulsion system
    :return: isp [s] Specific impulse of the propulsion system
    """

    return v_eff / g0


def veff_thrust_flow_rate(thrust: float, flow_rate: float) -> float:
    """
    Calculate the propulsion system's effective exit velocity

    :param thrust: [m*kg/s²] Thrust of the propulsion system
    :param flow_rate: [kg/s] Propellant flow rate of the propulsion system
    :return: veff [m/s] Effective exit velocity of the propulsion system
    """

    return thrust / flow_rate


def isp_thrust_flow_rate(g0: float, thrust: float, flow_rate: float) -> float:
    """
    Calculate the specific impulse of the propulsion system

    :param g0: [m/s²] Standard acceleration of gravity
    :param thrust: [m*kg/s²] Thrust of the propulsion system
    :param flow_rate: [kg/s] Propellant flow rate of the propulsion system
    :return: isp [s] Specific impulse of the propulsion system
    """

    return thrust / (g0 * flow_rate)


def flow_rate_isp_thrust(g0: float, isp: float, thrust: float) -> float:
    """
    Calculate the propellant flow rate of the propulsion system

    :param g0: [m/s²] Standard acceleration of gravity
    :param isp: [s] Specific impulse of the propulsion system
    :param thrust: [m*kg/s²] Thrust of the propulsion system
    :return: ṁ [kg/s] Propellant flow rate of the propulsion system
    """

    return thrust / (g0 * isp)


def thrust_isp_flow_rate(g0: float, isp: float, flow_rate: float) -> float:
    """
    Calculate the thrust of the propulsion system

    :param g0: [m/s²] Standard acceleration of gravity
    :param isp: [s] Specific impulse of the propulsion system
    :param flow_rate: [kg/s] Propellant flow rate of the propulsion system
    :return: T [m*kg/s²] Thrust of the propulsion system
    """

    return g0 * isp * flow_rate


# TODO: Calculating total thrust/effective velocity from multiple parallel propulsion systems
def veff_total(v_effs: List[float], flow_rates: List[float]) -> float:
    """
    Calculate the total effective exit velocity of a parallel propulsion system

    :param v_effs: *[m/s] List containing the effective exit velocities of the single propulsion systems
    :param flow_rates: *[kg/s] List containing the propellant flow rates of the single propulsion systems
    :return: v_eff [m/s] Total effective exit velocity of all parallel propulsion systems
    """

    assert len(v_effs) == len(flow_rates), "List sizes must be equal"

    v_effs_array = np.array(v_effs)
    flow_rates_array = np.array(flow_rates)

    return np.sum(v_effs_array * flow_rates_array) / np.sum(flow_rates_array)


def veff_veff_total(v_eff_total: float, v_effs: List[float], flow_rates: List[float]) -> float:
    """
    Calculate the effective exit velocity of one propulsion system in a set of parallel propulsion systems

    :param v_eff_total: [m/s] Total effective exit velocity of the parallel propulsion system
    :param v_effs: *[m/s] List containing the effective exit velocities of the single propulsion systems
        (except for the wanted value)
    :param flow_rates: *[kg/s] List containing the propellant flow rates of the single propulsion systems,
        with the first entry corresponding to the wanted effective exit velocity
    :return: v_eff [m/s] Effective exit velocity of the single propulsion system
    """

    assert len(flow_rates) - len(v_effs) == 1, \
        "There must be exactly one more entry in the flow_rates list than the v_effs list"

    v_effs_array = np.array(v_effs)
    flow_rates_array = np.array(flow_rates)

    return (v_eff_total * np.sum(flow_rates_array) - np.sum(v_effs_array * flow_rates_array[1:])) / flow_rates_array[0]


def flow_rate_veff_total(v_eff_total: float, v_effs: List[float], flow_rates: List[float]) -> float:
    """
    Calculate the propellant flow rate of one propulsion system in a set of parallel propulsion systems

    :param v_eff_total: [m/s] Total effective exit velocity of the parallel propulsion system
    :param v_effs: *[m/s] List containing the effective exit velocities of the single propulsion systems,
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
