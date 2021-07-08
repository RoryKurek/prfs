from . import Q_, ureg
import thermo as th
from .util import flash_to_VF


@ureg.check(None, '[pressure]', None, '[]', '[]', '[area]', None, '[]', None)
def api521_fire_wetted(flasher: th.flash.Flash,
                       P: Q_,
                       zs: list[float],
                       initial_VF: Q_,
                       final_VF: Q_,
                       A: Q_,
                       adequate_drainage: bool,
                       F: Q_ = Q_('1.0'),
                       air_cooler: bool = False) -> dict[str, Q_]:
    """
    Calculates the rate of vaporization for a vessel containing liquid
    under fire.

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

    The vaporization rate is then calculated based on
    :math:`\\dot{m} = Q \\times \\Delta H_{vap}`. API 521 does not
    specify exactly how the heat of vaporization is to be determined.
    The overall heat of vaporization of the liquid is not necessarily
    conservative - the "instantaneous" heat of vaporization may be
    higher. Therefore, a heat of vaporization is typically calculated
    over some interval (i.e. vaporization of the first 5 mol% of the
    liquid). Currently, this implementation assumes that the entire
    duty goes into vaporization, and not into the temperature increase
    of the bulk fluid. In other words:

    .. math:: \\Delta H_{vap} = (H(P, VF_{final}) - H(P, VF_{initial}))
              - C_{p, avg} (T_{final} - T_{initial})

    The average heat capacity is simply the average of the bulk heat
    capacities at the initial and final vapor fractions. This could
    potentially be improved in the future by integrating the heat
    capacity over the temperature interval - however, this would
    require a model of the bulk heat capacity vs temperature.

    .. note:: This function should typically be used for each equipment
        item individually. Attempting to use the function with the total
        wetted area of all equipment items will often lead to
        under-predicting the total duty (due to the wetted area exponent of
        0.82).

    Parameters
    ----------
    flasher : thermo.flash.Flash
        A flasher modeling the vessel liquid contents.
    P : pint.Quantity
        The relief pressure.
    zs : list[float]
        The mole fractions of the components.
    initial_VF, final_VF : float, float
        The vapor fractions bounding the interval over which the
        vaporization is to be considered.
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
    results : dict[str, pint.Quantity]
        Contains the results of the calculation - in particular, the
        vaporization rate `results['n']` and the heat duty
        `results['Q']`.

    """
    if final_VF <= initial_VF:
        raise ValueError('final_VF must be greater than initial_VF')

    if adequate_drainage:
        C = Q_('21000.0 BTU/hr/ft^2')
    else:
        C = Q_('34500.0 BTU/hr/ft^2')

    if air_cooler:
        E = 1.0
    else:
        E = 0.82

    # Conversion is to avoid unit strangeness due to [area] ^ 0.82
    Q = C * F * Q_((A**E).magnitude, 'ft^2')

    initial_state = flash_to_VF(flasher=flasher, P=P, VF=initial_VF, zs=zs)
    final_state = flash_to_VF(flasher=flasher, P=P, VF=final_VF, zs=zs)

    # TODO: Investigate integral of Cp
    avg_Cp = Q_((initial_state.Cp() + final_state.Cp()) / 2.0, 'J/K/mol')

    initial_T, final_T = initial_state.T, final_state.T

    # Calculate latent heat over the desired interval
    # TODO: Add option to exclude specific heat input or not
    interval_total_dH = Q_(final_state.H() - initial_state.H(), 'J/mol')
    interval_specific_dH = avg_Cp * Q_(final_T - initial_T, 'K')
    interval_latent_dH = interval_total_dH - interval_specific_dH
    latent_dH_per_vapor = interval_latent_dH / \
        (final_state.VF - initial_state.VF)

    # Calculate vapor generation based on heat input
    n = Q / latent_dH_per_vapor
    return dict(Q=Q,
                C=C,
                n=n,
                avg_Cp=avg_Cp,
                initial_T=initial_T,
                final_T=final_T,
                interval_total_dH=interval_total_dH,
                interval_specific_dH=interval_specific_dH,
                interval_latent_dH=interval_latent_dH,
                latent_dH_per_vapor=latent_dH_per_vapor)
