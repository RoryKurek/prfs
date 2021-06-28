from . import ureg, Q_
from math import pi


@ureg.check('[]', '[length]', '[length]', '[mass]*[length]**-3', None, None, None)
def darcy_dp(f, L, D, rho, v=None, Q=None, m=None) -> Q_:
    num_flow_specs = 0
    if v is not None:
        num_flow_specs += 1
    if Q is not None:
        v = Q / (pi / 4 * D**2)
        num_flow_specs += 1
    if m is not None:
        v = (m / rho) / (pi / 4 * D**2)
        num_flow_specs += 1

    if num_flow_specs > 1:
        raise Exception('Flow arguments are overspecified. Function requires exactly '
                        'one of the following: 1) v (velocity), 2) Q (volumetric '
                        'flow) or 3) m (mass flow)')
    if num_flow_specs < 1:
        raise Exception('Flow arguments are underspecified. Function requires exactly '
                        'one of the following: 1) v (velocity), 2) Q (volumetric '
                        'flow) or 3) m (mass flow)')

    return f * rho * L * v * v / 2 / D

