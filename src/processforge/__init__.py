"""
ProcessForge - A Python process simulation framework for chemical engineering.

Provides steady-state and dynamic process simulation capabilities including:
- Flowsheet modeling with recycle loop support
- Unit operations (Pump, Valve, Tank, Pipes, Strainer, Flash, Heater)
- Thermodynamic property calculations via CoolProp
- JSON-schema validated flowsheet configurations
- Results export to Zarr and Excel validation reports
- Equation-oriented (EO) steady-state solver via EOFlowsheet
"""

from .flowsheet import Flowsheet
from .solver import Solver
from .thermo import get_enthalpy_molar, get_Cp_molar, get_K_values, rachford_rice
from .validate import validate_flowsheet
from .result import (
    generate_validation_excel,
    plot_results,
    plot_timeseries,
    save_results_zarr,
)
from .units.pump import Pump
from .units.valve import Valve
from .units.strainer import Strainer
from .units.tank import Tank
from .units.pipes import Pipes
from .units.flash import Flash
from .units.heater import Heater
from .eo import EOFlowsheet, EOSolver

__version__ = "0.1.0"

__all__ = [
    "Flowsheet",
    "EOFlowsheet",
    "EOSolver",
    "Solver",
    "get_enthalpy_molar",
    "get_Cp_molar",
    "get_K_values",
    "rachford_rice",
    "validate_flowsheet",
    "generate_validation_excel",
    "plot_results",
    "plot_timeseries",
    "save_results_zarr",
    "Pump",
    "Valve",
    "Strainer",
    "Tank",
    "Pipes",
    "Flash",
    "Heater",
    "__version__",
]
