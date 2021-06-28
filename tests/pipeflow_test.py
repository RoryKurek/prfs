from psv_tools import Q_
from psv_tools.pipeflow import darcy_dp
import pytest


class TestDarcyDP:
    @pytest.mark.parametrize(
        'kwargs, dp', [
            ({'flow': Q_('5.0 ft/s'), 'f': Q_('0.016'), 'L': Q_('10.0 ft'),
              'D': Q_('2.0 in'), 'rho': Q_('62.4 lb/ft^3')}, Q_('0.161620941 psi')),
            ({'flow': Q_('0.10908307 ft^3/s'), 'f': Q_('0.016'), 'L': Q_('10.0 ft'),
              'D': Q_('2.0 in'), 'rho': Q_('62.4 lb/ft^3')}, Q_('0.161620941 psi')),
            ({'flow': Q_('24504.4226 lb/hr'), 'f': Q_('0.016'), 'L': Q_('10.0 ft'),
              'D': Q_('2.0 in'), 'rho': Q_('62.4 lb/ft^3')}, Q_('0.161620941 psi')),
        ]
    )
    def test_different_flow_inputs(self, kwargs, dp):
        assert darcy_dp(**kwargs).to('psi').magnitude == pytest.approx(dp.to('psi').magnitude)

    @pytest.mark.parametrize(
        'kwargs', [
            {'flow': Q_('5.0 mol/hr'), 'f': Q_('0.016'), 'L': Q_('10.0 ft'),
             'D': Q_('2.0 in'), 'rho': Q_('62.4 lb/ft^3')},
        ]
    )
    def test_fails_with_invalid_flow(self, kwargs):
        with pytest.raises(ValueError) as e:
            darcy_dp(**kwargs)
