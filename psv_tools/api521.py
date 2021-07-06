from collections import namedtuple
from . import ureg, Q_
import thermo as th
from scipy.optimize import root_scalar
from functools import cache


@ureg.wraps(None, (ureg.ft ** 2, None, ureg.dimensionless, None))
def fire_wetted_Q(A: Q_, adequate_drainage: bool, F: Q_ = Q_('1.0'),
                  air_cooler: bool = False) -> tuple[Q_, Q_]:
    """
    Calculates the heat duty for a vessel containing liquid under fire.

    The heat duty for non-air-cooler equipment items is calculated based on
    equations (7) and (8) from API Standard 521, 7th Ed., §4.4.13.2.4.2:

    .. math:: Q = C \\times F \\times A^{0.82}

    When the equipment in question is an air-cooled heat exchanger in
    liquid-cooling service, the wetted area exponent is changed to 1.0 per
    equations (21) and (22) in §4.4.13.2.8.4 of the same standard:

    .. math:: Q = C \\times F \\times A

    C is a constant dependent on whether the adequate drainage and
    firefighting is present. Per §4.4.13.2.4.2, "The determination of what
    constitutes adequate drainage is subjective and left to the user to
    decide but it should be designed to carry flammable/combustible liquids
    away from a vessel."

    .. note:: This function should typically be used for each equipment
        item individually. Attempting to use the function with the total
        wetted area of all equipment items will often lead to
        under-predicting the total duty (due to the wetted area exponent of
        0.82).

    Parameters
    ----------
    A : pint.Quantity
        The wetted area of the equipment.
    adequate_drainage : bool
        Whether or not the equipment has 'adequate' drainage and firefighting.
    F : pint.Quantity
        Environment factor to account for fire-proof insulation or
        earth-covered storage. Defaults to 1.0 for uninsulated equipment.
        See Table 5 in API Standard 521, 7th Ed., §4.4.13.2.4.2 for details.
    air_cooler : bool
        Whether or not the equipment is an air-cooled heat exchanger.

    Returns
    -------
    Q : pint.Quantity
        The calculated heat duty.
    C : pint.Quantity
        The duty constant used in the calculation.
    """
    if adequate_drainage:
        C = 21000.0
    else:
        C = 34500.0

    if air_cooler:
        E = 1.0
    else:
        E = 0.82

    return Q_(C * F * A ** E, 'BTU/hr'), Q_(C, 'BTU/hr/ft^2')


@ureg.wraps(None, (None, ureg.Pa, ureg.dimensionless, None))
@cache
def flash_to_VF(flasher: th.flash.Flash, P: Q_, VF: Q_, zs) \
        -> th.equilibrium.EquilibriumState:

    if VF == 0.0 or VF == 1.0:
        return flasher.flash(P=P, VF=VF, zs=zs)
    else:
        sat_liquid = flash_to_VF(flasher, P=Q_(P, 'Pa'), VF=Q_('0.0'), zs=zs)
        sat_vapor = flash_to_VF(flasher, P=Q_(P, 'Pa'), VF=Q_('1.0'), zs=zs)

        def check_T(T_val):
            return flasher.flash(T=T_val, P=P, zs=zs).VF - VF

        T = root_scalar(check_T, bracket=[sat_liquid.T, sat_vapor.T]).root
        return flasher.flash(T=T, P=P, zs=zs)


FireWettedResults = namedtuple('FireWettedResults', 'Q, C, n, avg_Cp, initial_T, final_T, '
                                                    'interval_total_dH, interval_specific_dH, '
                                                    'interval_latent_dH, latent_dH_per_vapor')


def fire_wetted(flasher: th.flash.Flash, P: Q_, zs, initial_VF: Q_,
                final_VF: Q_, A: Q_, adequate_drainage: bool, F: Q_ = Q_('1.0'),
                air_cooler: bool = False) -> FireWettedResults:
    """
    Calculates the heat duty for a vessel containing liquid under fire.

    The heat duty for non-air-cooler equipment items is calculated based on
    equations (7) and (8) from API Standard 521, 7th Ed., §4.4.13.2.4.2:

    .. math:: Q = C \\times F \\times A^{0.82}

    When the equipment in question is an air-cooled heat exchanger in
    liquid-cooling service, the wetted area exponent is changed to 1.0 per
    equations (21) and (22) in §4.4.13.2.8.4 of the same standard:

    .. math:: Q = C \\times F \\times A

    C is a constant dependent on whether the adequate drainage and
    firefighting is present. Per §4.4.13.2.4.2, "The determination of what
    constitutes adequate drainage is subjective and left to the user to
    decide but it should be designed to carry flammable/combustible liquids
    away from a vessel."

    .. note:: This function should typically be used for each equipment
        item individually. Attempting to use the function with the total
        wetted area of all equipment items will often lead to
        under-predicting the total duty (due to the wetted area exponent of
        0.82).

    Parameters
    ----------
    A : pint.Quantity
        The wetted area of the equipment.
    adequate_drainage : bool
        Whether or not the equipment has 'adequate' drainage and firefighting.
    F : pint.Quantity
        Environment factor to account for fire-proof insulation or
        earth-covered storage. Defaults to 1.0 for uninsulated equipment.
        See Table 5 in API Standard 521, 7th Ed., §4.4.13.2.4.2 for details.
    air_cooler : bool
        Whether or not the equipment is an air-cooled heat exchanger.

    Returns
    -------
    FireWettedResults

    """
    assert P.check('[pressure]') and initial_VF.check('[]') and \
           final_VF.check('[]') and A.check('[area]') and \
           final_VF > initial_VF

    if adequate_drainage:
        C = Q_('21000.0 BTU/hr/ft^2')
    else:
        C = Q_('34500.0 BTU/hr/ft^2')

    if air_cooler:
        E = 1.0
    else:
        E = 0.82

    # Conversion is to avoid unit strangeness due to [area] ^ 0.82
    Q = C * F * Q_((A ** E).magnitude, 'ft^2')

    initial_state = flash_to_VF(flasher=flasher, P=P, VF=initial_VF, zs=zs)
    final_state = flash_to_VF(flasher=flasher, P=P, VF=final_VF, zs=zs)

    # Currently uses average heat capacity based on the two endpoints
    # TODO: Investigate integral of dCp/dT
    avg_Cp = Q_((initial_state.Cp() + final_state.Cp()) / 2.0, 'J/K/mol')

    initial_T, final_T = initial_state.T, final_state.T

    # Calculate latent heat over the desired interval
    # TODO: Add option to exclude specific heat input or not
    interval_total_dH = Q_(final_state.H() - initial_state.H(), 'J/mol')
    interval_specific_dH = avg_Cp * Q_(final_T - initial_T, 'K')
    interval_latent_dH = interval_total_dH - interval_specific_dH
    latent_dH_per_vapor = interval_latent_dH / (final_state.VF - initial_state.VF)

    # Calculate vapor generation based on heat input
    n = Q / latent_dH_per_vapor
    return FireWettedResults(Q=Q, C=C, n=n, avg_Cp=avg_Cp, initial_T=initial_T,
                             final_T=final_T, interval_total_dH=interval_total_dH,
                             interval_specific_dH=interval_specific_dH,
                             interval_latent_dH=interval_latent_dH,
                             latent_dH_per_vapor=latent_dH_per_vapor)
