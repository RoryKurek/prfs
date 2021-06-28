import thermo as th
import thermo.interaction_parameters as ip


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
