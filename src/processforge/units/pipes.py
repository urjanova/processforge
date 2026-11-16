import copy

from loguru import logger

from ..eo.units.pipes_eo import PipesEOMixin
from .base import BaseUnitMixin
from .hydraulics import laminar_pipe_flow, resolve_properties
from .registry import register_unit


def _first_element(value):
    """Take a representative scalar from a snapshot-or-timeseries value."""
    if isinstance(value, list):
        return value[0] if value else None
    if isinstance(value, dict):
        return {k: _first_element(v) for k, v in value.items()}
    return value


class Pipes(BaseUnitMixin, PipesEOMixin):
    """
    Simple steady-state pipes model.
    Represents pipes connecting units, with a fixed pressure drop due to friction.
    Assumes isothermal (no temperature change).

    The auxiliary volumetric ``flow`` estimate uses the Hagen–Poiseuille
    relation with the mixture viscosity from the attached provider (or an
    explicit ``viscosity`` unit parameter [Pa·s]) and the ``length`` unit
    parameter (default 1 m). When no viscosity is available the estimate is
    set to 0 and a warning is logged — no hard-coded fluid default is used.
    """

    def __init__(self, name, delta_p=1000.0, diameter=0.1, length=1.0,
                 viscosity=None, **kwargs):
        self.name = name
        self.delta_p = delta_p
        self.diameter = diameter
        self.length = length
        self.viscosity = viscosity
        self.params = {
            "delta_p": delta_p,
            "diameter": diameter,
            "length": length,
            "viscosity": viscosity,
            **kwargs,
        }
        self._viscosity_warned = False

    def _mixture_viscosity(self, inlet):
        """Resolve the mixture viscosity [Pa·s] from provider or unit param."""
        provider = getattr(self, "_provider", None)
        props = resolve_properties(
            self.name,
            {
                "z": _first_element(inlet.get("z", {})),
                "T": _first_element(inlet.get("T")),
                "P": _first_element(inlet.get("P", 101325.0)),
            },
            {"viscosity": self.viscosity},
            getattr(provider, "get_thermo_properties", None) if provider is not None else None,
        )
        return props["viscosity"]

    def _run_impl(self, inlet):
        """
        Processes the inlet stream through a pipe unit, applying the fixed
        pressure drop (floored at 1000.0 Pa) and optionally estimating the
        volumetric flow via Hagen–Poiseuille. Temperature passes through
        unchanged.

        Handles both single-value snapshots and timeseries (lists).

        Args:
            inlet (dict): Inlet stream properties ('P' defaults to 101325.0 Pa).

        Returns:
            dict: Outlet stream with updated 'P', optional 'flow', and
            'unit' set to "Pipes".
        """
        outlet = copy.deepcopy(inlet)
        inlet_p = inlet.get("P", 101325.0)

        # Handle both scalar and list inputs
        if isinstance(inlet_p, list):
            outlet["P"] = [max(p - self.delta_p, 1000.0) for p in inlet_p]
            delta_p_used = [inlet_p[i] - outlet["P"][i] for i in range(len(inlet_p))]
        else:
            outlet["P"] = max(inlet_p - self.delta_p, 1000.0)
            delta_p_used = inlet_p - outlet["P"]

        viscosity = self._mixture_viscosity(inlet)
        if viscosity is not None:
            if isinstance(inlet_p, list):
                outlet["flow"] = [
                    laminar_pipe_flow(self.diameter, dp, viscosity, self.length)
                    for dp in delta_p_used
                ]
            else:
                outlet["flow"] = laminar_pipe_flow(
                    self.diameter, delta_p_used, viscosity, self.length
                )
        else:
            if not self._viscosity_warned:
                logger.warning(
                    f"Pipes '{self.name}': no mixture viscosity available from "
                    "the provider or unit parameters — setting the auxiliary "
                    "'flow' estimate to 0. Add \"viscosity\" to the unit or use "
                    "a provider that supplies it."
                )
                self._viscosity_warned = True
            outlet["flow"] = [0.0] * len(inlet_p) if isinstance(inlet_p, list) else 0.0

        outlet["T"] = inlet["T"]
        outlet["unit"] = "Pipes"
        return outlet

    def run_dynamic(self, inlet_stream_ts, time_range, t_eval, solver):
        """
        Simulates dynamic operation of the pipes unit over a time step dt.
        For this simple model, assumes no accumulation and behaves like steady-state.

        Parameters:
            inlet_stream_ts: Inlet stream timeseries dictionary.
            time_range: Time range for simulation (not used in this simple model).
            t_eval: Time points for evaluation (not used in this simple model).
            solver: Solver object (not used in this simple model).

        Returns:
            dict: Outlet stream dictionary, same as run method.
        """
        return self.run(inlet_stream_ts)


register_unit("Pipes", Pipes)
