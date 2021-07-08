from scipy.optimize import root_scalar
from functools import cache
import thermo as th
import thermo.interaction_parameters as ip
from . import ureg, Q_


def create_VL_flasher(names: list[str]) -> th.flash.FlashVL:
    # Look up chemical constants/properties from database
    constants, properties = th.ChemicalConstantsPackage.from_IDs(names)

    # Look up binary interaction parameters from database
    kijs = ip.IPDB.get_ip_asymmetric_matrix('ChemSep PR', constants.CASs, 'kij')

    eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs,
                  'omegas': constants.omegas, 'kijs': kijs}

    # Create individual gas and liquid phase EOSs
    gas = th.CEOSGas(th.PRMIX, eos_kwargs,
                     HeatCapacityGases=properties.HeatCapacityGases)
    liquid = th.CEOSLiquid(th.PRMIX, eos_kwargs,
                           HeatCapacityGases=properties.HeatCapacityGases)

    # Create VL flasher
    return th.FlashVL(constants, properties, liquid=liquid, gas=gas)


@ureg.wraps(None, (None, ureg.Pa, ureg.dimensionless, None))
@cache
def flash_to_VF(flasher: th.flash.Flash, P: Q_, VF: Q_, zs) \
        -> th.equilibrium.EquilibriumState:

    if VF == 0.0 or VF == 1.0:
        return flasher.flash(P=P, VF=VF, zs=zs)
    else:
        sat_liquid = flash_to_VF(flasher, P=Q_(P, 'Pa'), VF=Q_('0.0'), zs=zs)
        sat_vapor = flash_to_VF(flasher, P=Q_(P, 'Pa'), VF=Q_('1.0'), zs=zs)

        def check_T(T_val):
            return flasher.flash(T=T_val, P=P, zs=zs).VF - VF

        T = root_scalar(check_T, bracket=[sat_liquid.T, sat_vapor.T]).root
        return flasher.flash(T=T, P=P, zs=zs)
