from . import ureg, Q_
import thermo as th
from scipy.optimize import root_scalar


@ureg.wraps(None, (ureg.ft ** 2, None, ureg.dimensionless, None))
def fire_wetted_Q(A: Q_, adequate_drainage: bool, F: Q_ = Q_('1.0'),
                  air_cooler: bool = False) -> dict[str, Q_]:
    """Calculates the heat absorbed by equipment containing liquid under fire
    based on API 521.

    The heat duty is normally calculated based on equations (7) and (8) from
    API Standard 521, 7th Ed., §4.4.13.2.4.2:

    .. math:: Q = C \times F \times A ^ 0.82

    When the equipment in question is an air-cooled heat exchanger in
    liquid-cooling service, the wetted area exponent is changed to 1.0 per
    equations (21) and (22) in §4.4.13.2.8.4 of the same standard:

    .. math:: Q = C \times F \times A

    C is a constant dependent on whether the adequate drainage and
    firefighting is present. Per §4.4.13.2.4.2, "The determination of what
    constitutes adequate drainage is subjective and left to the user to decide
    but it should be designed to carry flammable/combustible liquids away from
    a vessel."

    Environment factor F is used to account for fire-proof insulation or
    earth-covered storage. See Table 5 in §4.4.13.2.4.2 for details.

    :param A: Wetted surface area of the equipment item as a pint quantity
    with [area] dimensionality
    :param adequate_drainage: Whether or not the equipment
    :param F:
    :param air_cooler:
    :return:
    """

    if adequate_drainage:
        C = 21000.0
    else:
        C = 34500.0

    if air_cooler:
        E = 1.0
    else:
        E = 0.82

    return {'C': Q_(C, 'BTU/hr/ft^2'),
            'Q': Q_(C * F * A ** E, 'BTU/hr')}


@ureg.wraps(None, (None, ureg.Pa, ureg.dimensionless, ureg.dimensionless, None))
def flash_to_VFs(flasher: th.flash.Flash, P: Q_, initial_VF: Q_, final_VF: Q_,
                zs) -> tuple[th.equilibrium.EquilibriumState, th.equilibrium.EquilibriumState]:
    assert initial_VF != 1.0 and final_VF != 0.0 and initial_VF < final_VF

    sat_liquid = flasher.flash(P=P, VF=0, zs=zs)
    sat_vapor = flasher.flash(P=P, VF=1, zs=zs)

    def check_T(T, VF):
        return flasher.flash(T=T, P=P, zs=zs).VF - VF

    if initial_VF == 0.0:
        initial_state = sat_liquid
    else:
        T = root_scalar(check_T, args=initial_VF, bracket=[sat_liquid.T, sat_vapor.T]).root
        initial_state = flasher.flash(T=T, P=P, zs=zs)

    if final_VF == 1.0:
        final_state = sat_vapor
    else:
        T = root_scalar(check_T, args=final_VF, bracket=[sat_liquid.T, sat_vapor.T]).root
        final_state = flasher.flash(T=T, P=P, zs=zs)

    return initial_state, final_state


def fire_wetted(flasher: th.flash.Flash, P: Q_, zs, initial_VF: Q_,
                final_VF: Q_, A: Q_, adequate_drainage: bool, F: Q_ = Q_('1.0'),
                air_cooler: bool = False):
    assert P.check('[pressure]') and initial_VF.check('[]') and \
           final_VF.check('[]') and A.check('[area]')

    results = fire_wetted_Q(A=A, adequate_drainage=adequate_drainage, F=F,
                            air_cooler=air_cooler)

    initial_state, final_state = \
        flash_to_VFs(flasher=flasher, P=P, initial_VF=initial_VF, final_VF=final_VF, zs=zs)

    results['avg_Cp'] = Q_((initial_state.Cp() + final_state.Cp()) / 2.0, 'J/K/mol')
    # TODO: Investigate integral of dCp/dT

    results['interval_total_dH'] = Q_(final_state.H() - initial_state.H(), 'J/mol')
    results['interval_specific_dH'] = results['avg_Cp'] * Q_(final_state.T - initial_state.T, 'K')
    results['interval_latent_dH'] = results['interval_total_dH'] - results['interval_specific_dH']
    results['latent_dH_per_vapor'] = results['interval_latent_dH'] / (final_state.VF - initial_state.VF)

    results['n'] = results['Q'] / results['latent_dH_per_vapor']
    return results
