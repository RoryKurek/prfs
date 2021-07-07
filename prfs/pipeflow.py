from . import ureg, Q_
from math import pi


@ureg.check(None, '[]', '[length]', '[length]', '[mass]/[volume]')
def darcy_dp(flow: Q_, f: Q_, L: Q_, D: Q_, rho: Q_) -> Q_:
    if flow.check('[length]/[time]'):
        v = flow
    elif flow.check('[volume]/[time]'):
        v = flow / (pi * D ** 2 / 4)
    elif flow.check('[mass]/[time]'):
        v = flow / rho / (pi * D ** 2 / 4)
    else:
        raise ValueError('flow argument must be a mass flow, volumetric flow '
                         'or velocity.')

    return f * rho * L * v * v / 2 / D

