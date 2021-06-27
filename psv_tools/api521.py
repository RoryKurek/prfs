from . import ureg, Q_


@ureg.wraps(None, (ureg.ft ** 2, None, ureg.dimensionless))
def fire_wetted_Q(A: Q_, adequate_drainage: bool, F = 1.0) -> Q_:
    if adequate_drainage:
        C = 21000.0
    else:
        C = 34500.0

    return {'C': Q_(C, 'BTU/hr/ft^2'),
            'Q': Q_(C * F * A ** 0.82, 'BTU/hr')}
