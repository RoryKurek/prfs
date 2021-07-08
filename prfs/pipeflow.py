from typing import Optional
from . import ureg, Q_
from math import pi


@ureg.check('[length]/[time]', '[mass]/[volume]', '[]', '[length]', '[length]', '[]', '[length]')
def segment_dP_incompressible(v: Q_, rho_mass: Q_, f: Q_, D: Q_, L: Q_ = Q_('0.0 ft'), K: Q_ = Q_('0.0'), dz: Q_ = Q_('0.0 ft')) -> Q_:
    """
    Calculates the pressure drop across a pipe segment of constant
    diameter for a given flow velocity. Assumes incompressible flow;
    therefore, the velocity and density are assumed to be assumed to be
    constant between inlet and outlet.

    Flow is calculated using Bernoulli's equation, with head loss based
    on the Weisbach formula.

    .. math:: \\Delta P = \\rho \\left[ \\left( f \\frac{L}{D} + K \\right)
              \\frac{v^2}{2D} + \\Delta z \\right]

    Parameters
    ----------
    v
        Fluid velocity (assumed constant between inlet and outlet
    rho_mass
        Fluid mass density (assumed constant between inlet and outlet
    f
        Fluid Darcy friction factor
    D
        Piping diameter (assumed constant)
    L
        Length of piping
    K
        Head loss coefficient for minor losses (i.e. fittings)
    dz
        Elevation change between inlet and outlet

    Returns
    -------
    dP
        Pressure drop between the inlet and outlet
    """
    return rho_mass * ((f * L / D + K) * v * v / 2 + dz)


@ureg.check(None, '[length]', None)
def calculate_velocity(flow: Q_, D: Q_, rho_mass: Optional[Q_] = None):
    """
    Calculate the average velocity of flow through a pipe given a mass
    flow or volumetric flow and the flow diameter. If the specified
    flow is a mass flow rate, a mass density is also required.

    If ``flow`` is a velocity, it is immediately returned.

    .. math:: v &= \\frac{\\dot{V}}{\\pi D^2 / 4} \\\\
                &= \\frac{\\dot{m} / \\rho}{\\pi D^2 / 4}

    Parameters
    ----------
    flow
        Flow as a velocity, mass flow rate or volumetric flow rate
    D
        Flow diameter
    rho_mass
        Fluid mass density

    Returns
    -------
    v
        Average fluid velocity
    """
    if flow.check('[length]/[time]'):
        return flow
    elif flow.check('[volume]/[time]'):
        return flow / (pi * D ** 2 / 4)
    elif flow.check('[mass]/[time]'):
        if rho_mass is not None and rho_mass.check('[mass]/[volume]'):
            return flow / rho_mass / (pi * D ** 2 / 4)
        else:
            raise ValueError('If `flow` is a mass flow rate, optional argument '
                             '`rho_mass` must be passed and must be a mass density')
    # TODO: Add handling for molar flows
    else:
        raise ValueError('flow argument must be a mass flow, volumetric flow '
                         'or velocity.')
