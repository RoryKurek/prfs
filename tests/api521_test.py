from psv_tools import Q_
from psv_tools.api521 import fire_wetted_Q
import pytest


class TestFireWettedQ:
    @pytest.mark.parametrize(
        'A, adequate_drainage, F, C, Q',
        [(Q_('1000.0 ft^2'), False, Q_('1.0'), Q_('34500.0 BTU/hr/ft^2'), Q_('9949908.6857 BTU/hr')),
         (Q_('500.0 m^2'), True, Q_('0.75'), Q_('21000.0 BTU/hr/ft^2'), Q_('5292056.03 W'))]
    )
    def test_works_with_quantity_inputs(self, A, adequate_drainage, F, C, Q):
        results = fire_wetted_Q(A, adequate_drainage, F)
        assert results['C'] == C
        assert results['Q'].to('BTU/hr').magnitude == pytest.approx(Q.to('BTU/hr').magnitude)

    def test_fails_with_non_quantity_area(self):
        with pytest.raises(ValueError):
            fire_wetted_Q(A=3500.0, adequate_drainage=True)
