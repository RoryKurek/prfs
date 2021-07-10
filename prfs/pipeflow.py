from math import pi


def segment_dP_incompressible(w: float, rho_mass: float, f: float, D: float, L: float = 0.0, K: float = 0.0,
                              dz: float = 0.0) -> float:
    """
    Calculates the pressure drop across a pipe segment of constant
    diameter for a given flow velocity. Assumes incompressible flow.

    Flow is calculated using Bernoulli's equation. As the flow is
    incompressible, the velocity and density are assumed to be assumed
    to be constant between inlet and outlet. The head loss component is
    based on the Weisbach formula.

    .. math:: \\Delta P = \\rho_m \\left[ \\left( f \\frac{L}{D} + K \\right)
              \\frac{w^2}{2} + g \\Delta z \\right]

    The gravitational acceleration :math:`g` is considered to be
    constant at 9.80665 m/s².

    Parameters
    ----------
    w
        Fluid velocity (assumed constant between inlet and outlet)
        [m/s]
    rho_mass
        Fluid mass density (assumed constant between inlet and outlet)
        [kg/m³]
    f
        Fluid Darcy friction factor [dimensionless]
    D
        Piping diameter (assumed constant) [m]
    L
        Length of piping [m]
    K
        Head loss coefficient for minor losses (i.e. fittings)
        [dimensionless]
    dz
        Elevation change between inlet and outlet [m]

    Returns
    -------
    dP
        Pressure drop between the inlet and outlet [Pa]
    """
    return rho_mass * ((f * L / D + K) * w * w / 2 + 9.80665 * dz)


def pipe_velocity_from_mass_flow(m_dot: float, D: float, rho_mass: float):
    """
    Calculate the average velocity of flow through a pipe of constant
    cross-section based on a mass flow rate.

    .. math:: w = \\frac{\\dot{m} / \\rho_m}{\\pi D^2 / 4}

    Parameters
    ----------
    m_dot
        Mass flow rate [kg/s]
    D
        Flow diameter [m]
    rho_mass
        Fluid mass density [kg/m³]

    Returns
    -------
    w
        Average fluid velocity [m/s]
    """
    return m_dot / rho_mass / (pi * D ** 2 / 4)
