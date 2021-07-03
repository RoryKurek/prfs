from psv_tools import Q_
from psv_tools.api521 import fire_wetted_Q
import pytest
from pint import DimensionalityError


class TestFireWettedQ:
    @pytest.mark.parametrize(
        'A, adequate_drainage, F, air_cooler, C, Q',
        [(Q_('1000.0 ft^2'), False, Q_('1.0'), False, Q_('34500.0 BTU/hr/ft^2'), Q_('9949908.6857 BTU/hr')),
         (Q_('500.0 m^2'), True, Q_('0.75'), False, Q_('21000.0 BTU/hr/ft^2'), Q_('5292056.03 W')),
         (Q_('500.0 ft^2'), False, Q_('1.0'), True, Q_('34500.0 BTU/hr/ft^2'), Q_('17250000.0 BTU/hr'))]
    )
    def test_works_with_quantity_inputs(self, A, adequate_drainage, F, C, air_cooler, Q):
        Q, C = fire_wetted_Q(A, adequate_drainage, F, air_cooler)
        assert C == C
        assert Q.to('BTU/hr').magnitude == pytest.approx(Q.to('BTU/hr').magnitude)

    @pytest.mark.parametrize(
        'A, adequate_drainage, F, air_cooler',
        [(1000.0, False, Q_('1.0'), False),
         (Q_('500.0 m^2'), True, 0.75, False)]
    )
    def test_fails_with_non_quantity_inputs(self, A, adequate_drainage, F, air_cooler):
        with pytest.raises(ValueError):
            fire_wetted_Q(A, adequate_drainage, F, air_cooler)

    @pytest.mark.parametrize(
        'A, adequate_drainage, F, air_cooler',
        [(Q_('500.0 lb'), True, Q_('0.75'), False),
         (Q_('1000.0 ft^2'), False, Q_('1.0 psi'), True)]
    )
    def test_fails_with_incorrect_units_on_inputs(self, A, adequate_drainage, F, air_cooler):
        with pytest.raises(DimensionalityError):
            fire_wetted_Q(A, adequate_drainage, F, air_cooler)