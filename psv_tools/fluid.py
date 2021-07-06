import thermo as th
from . import ureg
import inspect


METHOD_UNITS = {
    'A': ('J/mol', ()),
    'A_dep': ('J/mol', ()),
    'A_formation_ideal_gas': ('J/mol', (None,)),
    'A_ideal_gas': ('J/mol', (None,)),
    'A_mass': ('J/kg', (None,)),
    'A_reactive': ('J/mol', ()),
    'Bvirial': ('m^3/mol', (None,)),
    'Cp': ('J/mol/K', ()),
    'Cp_dep': ('J/mol/K', (None,)),
    'Cp_ideal_gas': ('J/mol/K', (None,)),
    'Cp_mass': ('J/kg/K', (None,)),
    'Cv': ('J/mol/K', ()),
    'Cv_dep': ('J/mol/K', (None,)),
    'Cv_ideal_gas': ('J/mol/K', (None,)),
    'Cv_mass': ('J/kg/K', (None,)),
    'G': ('J/mol', ()),
    'G_dep': ('J/mol', ()),
    'G_formation_ideal_gas': ('J/mol', (None,)),
    'G_ideal_gas': ('J/mol', (None,)),
    'G_mass': ('J/kg', (None,)),
    'G_reactive': ('J/mol', ()),
    'H': ('J/mol', ()),
    'H_dep': ('J/mol', (None,)),
    'H_formation_ideal_gas': ('J/mol', (None,)),
    'H_ideal_gas': ('J/mol', (None,)),
    'H_mass': ('J/kg', (None,)),
    'H_reactive': ('J/mol', ()),
    'Hc': ('J/mol', (None,)),
    'Hc_lower': ('J/mol', (None,)),
    'Hc_lower_mass': ('J/kg', (None,)),
    'Hc_lower_normal': ('J/m^3', (None,)),
    'Hc_lower_standard': ('J/m^3', (None,)),
    'Hc_mass': ('J/kg', (None,)),
    'Hc_normal': ('J/m^3', (None,)),
    'Hc_standard': ('J/m^3', (None,)),
    'Joule_Thomson': ('K/Pa', ()),
    'MW': ('g/mol', (None,)),
    'MWs': ('g/mol', (None,)),
    'Pmc': ('Pa', (None,)),
    'S': ('J/mol/K', ()),
    'S_dep': ('J/mol/K', (None,)),
    'S_formation_ideal_gas': ('J/mol/K', (None,)),
    'S_ideal_gas': ('J/mol/K', (None,)),
    'S_mass': ('J/kg/K', (None,)),
    'S_reactive': ('J/mol/K', ()),
    'Tmc': ('K', (None,)),
    'U': ('J/mol', ()),
    'U_dep': ('J/mol', ()),
    'U_formation_ideal_gas': ('J/mol', (None,)),
    'U_ideal_gas': ('J/mol', (None,)),
    'U_mass': ('J/kg', (None,)),
    'V': ('m^3/mol', ()),
    'V_dep': ('m^3/mol', ()),
    'V_gas': ('m^3/mol', (None,)),
    'V_gas_normal': ('m^3/mol', (None,)),
    'V_gas_standard': ('m^3/mol', (None,)),
    'V_ideal_gas': ('m^3/mol', (None,)),
    'V_iter': ('m^3/mol', (None, None)),
    'V_liquid_ref': ('m^3/mol', (None,)),
    'V_liquids_ref': ('m^3/mol', ()),
    'V_mass': ('m^3/kg', (None,)),
    'Vmc': ('m^3/mol', (None,)),
    'Wobbe_index': ('J/mol', (None,)),
    'Wobbe_index_lower': ('J/mol', (None,)),
    'Wobbe_index_lower_mass': ('J/kg', (None,)),
    'Wobbe_index_lower_normal': ('J/m^3', (None,)),
    'Wobbe_index_lower_standard': ('J/m^3', (None,)),
    'Wobbe_index_mass': ('J/kg', (None,)),
    'Wobbe_index_normal': ('J/m^3', (None,)),
    'Wobbe_index_standard': ('J/m^3', (None,)),
    'alpha': ('m^2/s', (None,)),
    'd2P_dT2': ('Pa/K^2', ()),
    'd2P_dT2_frozen': ('Pa/K^2', ()),
    'd2P_dTdV': ('mol*Pa^2/J/K', ()),
    'd2P_dTdV_frozen': ('mol*Pa^2/J/K', ()),
    'd2P_dV2': ('Pa*mol^2/m^6', ()),
    'd2P_dV2_frozen': ('Pa*mol^2/m^6', ()),
    'dA_dP': ('J/mol/Pa', ()),
    'dA_dP_T': ('J/mol/Pa', ()),
    'dA_dP_V': ('J/mol/Pa', ()),
    'dA_dT': ('J/mol/K', ()),
    'dA_dT_P': ('J/mol/K', ()),
    'dA_dT_V': ('J/mol/K', ()),
    'dA_dV_P': ('J/m^3', ()),
    'dA_dV_T': ('J/m^3', ()),
    'nu': ('m^2/s', (None,)),
    'pseudo_Pc': ('Pa', (None,)),
    'pseudo_Tc': ('K', (None,)),
    'pseudo_Vc': ('m^3/mol', (None,)),
    'rho_mass': ('kg/m^3', (None,)),
    'rho_mass_liquid_ref': ('kg/m^3', (None,))
}


PROPERTY_UNITS = {
    'Gfgs': 'J/mol',
    'Gfgs_mass': 'J/kg',
    'Hcs': 'J/mol',
    'Hcs_lower': 'J/mol',
    'Hcs_lower_mass': 'J/kg',
    'Hcs_mass': 'J/kg',
    'Hf_STPs': 'J/mol',
    'Hf_STPs_mass': 'J/kg',
    'Hfgs': 'J/mol',
    'Hfgs_mass': 'J/kg',
    'Hfus_Tms': 'J/mol',
    'Hfus_Tms_mass': 'J/kg',
    'Hsub_Tts': 'J/mol',
    'Hsub_Tts_mass': 'J/kg',
    'Hvap_298s': 'J/mol',
    'Hvap_298s_mass': 'J/kg',
    'Hvap_Tbs': 'J/mol',
    'Hvap_Tbs_mass': 'J/kg',
    'Pcs': 'Pa',
    'Psat_298s': 'Pa',
    'Pts': 'Pa',
    'RI_Ts': 'K',
    'S0gs': 'J/mol/K',
    'S0gs_mass': 'J/kg/K',
    'Sfgs': 'J/mol/K',
    'Sfgs_mass': 'J/kg/K',
    'Stockmayers': 'K',
    'Tautoignitions': 'K',
    'Tbs': 'K',
    'Tcs': 'K',
    'Tflashs': 'K',
    'Tms': 'K',
    'Tts': 'K',
    'Van_der_Waals_areas': 'm^2/mol',
    'Van_der_Waals_volumes': 'm^3/mol',
    'Vcs': 'm^3/mol',
    'Vmg_STPs': 'm^3/mol',
    'Vml_60Fs': 'm^3/mol',
    'Vml_STPs': 'm^3/mol',
    'Vml_Tms': 'm^3/mol',
    'Vms_Tms': 'm^3/mol',
    'conductivities': 'S/m',
    'conductivity_Ts': 'K',
}


def create_property(name, units=None):
    @ureg.wraps(units, None)
    def getter(self):
        return getattr(self._state, name)

    return property(getter)


def link_properties(cls):
    # Create properties in the StateUnitsWrapper class for all properties in EquilibriumState
    for name, member in inspect.getmembers(th.equilibrium.EquilibriumState):
        if not inspect.isfunction(member) and name[0] != '_':
            units = PROPERTY_UNITS.get(name, None)
            setattr(cls, name, create_property(name, units))

    return cls


@link_properties
class StateUnitsWrapper:
    def __init__(self, state: th.equilibrium.EquilibriumState):
        self._state = state

        for name, member in inspect.getmembers(state):
            if inspect.ismethod(member):
                if name in METHOD_UNITS:
                    setattr(self, name, ureg.wraps(*METHOD_UNITS[name])(member))
                else:
                    setattr(self, name, member)
