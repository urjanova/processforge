"""
ProcessForge - A Python process simulation framework for chemical engineering.

Provides steady-state and dynamic process simulation capabilities including:
- Flowsheet modeling with recycle loop support
- Unit operations (Pump, Valve, Tank, Pipes, Strainer, Flash, Heater)
- Thermodynamic property calculations via CoolProp
- JSON-schema validated flowsheet configurations
- Results export to CSV, JSON, and Excel validation reports
"""

from .flowsheet import Flowsheet
from .solver import Solver
from .thermo import get_enthalpy_molar, get_Cp_molar, get_K_values, rachford_rice
from .validate import validate_flowsheet
from .result import (
    save_results_csv,
    save_timeseries_csv,
    save_results_json,
    save_timeseries_json,
    generate_validation_excel,
)
from .units.pump import Pump
from .units.valve import Valve
from .units.strainer import Strainer
from .units.tank import Tank
from .units.pipes import Pipes
from .units.flash import Flash
from .units.heater import Heater

__version__ = "0.1.0"

__all__ = [
    "Flowsheet",
    "Solver",
    "get_enthalpy_molar",
    "get_Cp_molar",
    "get_K_values",
    "rachford_rice",
    "validate_flowsheet",
    "save_results_csv",
    "save_timeseries_csv",
    "save_results_json",
    "save_timeseries_json",
    "generate_validation_excel",
    "Pump",
    "Valve",
    "Strainer",
    "Tank",
    "Pipes",
    "Flash",
    "Heater",
    "__version__",
]
