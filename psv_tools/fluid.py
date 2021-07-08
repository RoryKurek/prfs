from itertools import repeat
import thermo as th
from . import ureg
import inspect, re
from pint import UndefinedUnitError


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

find_return_units = re.compile(r'Returns.+\[([^-]*)\]', re.DOTALL)



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

        for name, member in inspect.getmembers(state, inspect.ismethod):
            if name[0] != '_' and hasattr(member, '__doc__') and member.__doc__ is not None:
                match = find_return_units.search(member.__doc__)
                if match is None:
                    setattr(self, name, member)
                else:
                    return_units = match.group(1)
                    num_args = len(inspect.signature(member).parameters)
                    args = tuple(repeat(None, num_args))
                    try:
                        setattr(self, name, ureg.wraps(return_units, args)(member))
                    except UndefinedUnitError:
                        pass
