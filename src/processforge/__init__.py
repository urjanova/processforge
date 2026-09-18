"""
Processforge - A Python process simulation framework for chemical engineering.

Provides steady-state and dynamic process simulation capabilities including:
- Flowsheet modeling with recycle loop support
- Unit operations (Pump, Valve, Tank, Pipes, Strainer, Flash, Heater)
- Thermodynamic property calculations via CoolProp
- JSON-schema validated flowsheet configurations
- Results export to Zarr and Excel validation reports
- Equation-oriented (EO) steady-state solver via EOFlowsheet
"""

__version__ = "0.6.3"

from .eo import EOFlowsheet, EOSolver
from .flowsheet import Flowsheet
from .provenance import build_run_info
from .result import (
    save_results_zarr,
)
from .runner import (
    ApplyResult,
    ConvergenceError,
    FlowsheetValidationError,
    ProcessforgeRunError,
    ProviderRunError,
    ProviderUnavailableError,
    RunResult,
    StatePersistenceError,
    apply_flowsheet,
    run_flowsheet,
)
from .solver import Solver
from .thermo import rachford_rice
from .units.flash import Flash
from .units.heater import Heater
from .units.pipes import Pipes
from .units.pump import Pump
from .units.strainer import Strainer
from .units.tank import Tank
from .units.valve import Valve
from .utils.validate_flowsheet import validate_flowsheet

__all__ = [
    "ApplyResult",
    "ConvergenceError",
    "EOFlowsheet",
    "EOSolver",
    "Flash",
    "Flowsheet",
    "FlowsheetValidationError",
    "Heater",
    "Pipes",
    "ProcessforgeRunError",
    "ProviderRunError",
    "ProviderUnavailableError",
    "Pump",
    "RunResult",
    "Solver",
    "StatePersistenceError",
    "Strainer",
    "Tank",
    "Valve",
    "__version__",
    "apply_flowsheet",
    "build_run_info",
    "rachford_rice",
    "run_flowsheet",
    "save_results_zarr",
    "validate_flowsheet",
]
