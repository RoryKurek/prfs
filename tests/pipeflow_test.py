from psv_tools import Q_
from psv_tools.pipeflow import darcy_dp
import pytest


class TestDarcyDP:
    @pytest.mark.parametrize(
        'kwargs, dp', [
            ({'f': Q_('0.016'), 'L': Q_('10.0 ft'), 'D': Q_('2.0 in'),
              'rho': Q_('62.4 lb/ft^3'), 'v': Q_('5.0 ft/s')}, Q_('0.161620941 psi')),
            ({'f': Q_('0.016'), 'L': Q_('10.0 ft'), 'D': Q_('2.0 in'),
              'rho': Q_('62.4 lb/ft^3'), 'Q': Q_('0.10908307 ft^3/s')}, Q_('0.161620941 psi')),
            ({'f': Q_('0.016'), 'L': Q_('10.0 ft'), 'D': Q_('2.0 in'),
              'rho': Q_('62.4 lb/ft^3'), 'm': Q_('24504.4226 lb/hr')}, Q_('0.161620941 psi')),
        ]
    )
    def test_different_flow_inputs(self, kwargs, dp):
        assert darcy_dp(**kwargs).to('psi').magnitude == pytest.approx(dp.to('psi').magnitude)

    @pytest.mark.parametrize(
        'kwargs', [
            {'f': Q_('0.016'), 'L': Q_('10.0 ft'), 'D': Q_('2.0 in'),
             'rho': Q_('62.4 lb/ft^3'), 'v': Q_('5.0 ft/s'),
             'Q': Q_('0.10908307 ft^3/s')},
            {'f': Q_('0.016'), 'L': Q_('10.0 ft'), 'D': Q_('2.0 in'),
             'rho': Q_('62.4 lb/ft^3'), 'Q': Q_('0.10908307 ft^3/s'),
             'm': Q_('24504.4226 lb/hr')},
            {'f': Q_('0.016'), 'L': Q_('10.0 ft'), 'D': Q_('2.0 in'),
             'rho': Q_('62.4 lb/ft^3'), 'Q': Q_('0.10908307 ft^3/s'),
             'v': Q_('5.0 ft/s'), 'm': Q_('24504.4226 lb/hr')},
        ]
    )
    def test_inputs_overspecified(self, kwargs):
        with pytest.raises(Exception) as e:
            darcy_dp(**kwargs)
        assert str(e.value) == 'Flow arguments are overspecified. Function requires ' \
                               'exactly one of the following: 1) v (velocity), ' \
                               '2) Q (volumetric flow) or 3) m (mass flow)'

    @pytest.mark.parametrize(
        'kwargs', [
            {'f': Q_('0.016'), 'L': Q_('10.0 ft'), 'D': Q_('2.0 in'),
             'rho': Q_('62.4 lb/ft^3')},
        ]
    )
    def test_inputs_underspecified(self, kwargs):
        with pytest.raises(Exception) as e:
            darcy_dp(**kwargs)
        assert str(e.value) == 'Flow arguments are underspecified. Function requires ' \
                               'exactly one of the following: 1) v (velocity), ' \
                               '2) Q (volumetric flow) or 3) m (mass flow)'
