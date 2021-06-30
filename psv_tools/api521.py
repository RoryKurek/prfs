from . import ureg, Q_
import thermo as th
from scipy.optimize import root_scalar


@ureg.wraps(None, (ureg.ft ** 2, None, ureg.dimensionless, None))
def fire_wetted_Q(A: Q_, adequate_drainage: bool, F = 1.0,
                  air_cooler: bool = False) -> Q_:
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

# @ureg.wraps(None, (None, ureg.ft ** 2, None, ureg.dimensionless))
# def fire_wetted(flasher: th.flash.Flash, P_rel, VF, A, adequate_drainage: bool, F = 1.0):
#
