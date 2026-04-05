"""EOFlowsheet: equation-oriented steady-state flowsheet solver."""
from __future__ import annotations

from copy import deepcopy

import numpy as np
from loguru import logger

from .jacobian import GlobalJacobianManager
from .solver import EOSolver


def _get_inlets(unit_cfg: dict) -> list[str]:
    """Return inlet stream name(s) as a list."""
    inlet = unit_cfg.get("in")
    if isinstance(inlet, list):
        return inlet
    return [inlet] if inlet else []


def _get_outlets(unit_cfg: dict) -> list[str]:
    """Return outlet stream name(s) — handles Flash (out_vap + out_liq)."""
    if unit_cfg.get("type") == "Flash":
        outlets = []
        if unit_cfg.get("out_liq"):
            outlets.append(unit_cfg["out_liq"])
        if unit_cfg.get("out_vap"):
            outlets.append(unit_cfg["out_vap"])
        return outlets
    out = unit_cfg.get("out")
    return [out] if out else []


class EOFlowsheet:
    """
    Equation-oriented steady-state flowsheet solver.

    Reads the same JSON config as ``Flowsheet``.  Assembles all unit equations
    and stream connections into a global ``F(x) = 0`` system and solves
    simultaneously via Newton-Raphson.

    Recycles are handled natively — no tear stream initialisation needed.

    Args:
        config: Validated flowsheet configuration dict.
        backend: EO solver backend — ``"pyomo"`` (default), ``"casadi"``,
                 or ``"scipy"``.

    Raises:
        NotImplementedError: If the flowsheet contains a Tank unit (dynamic
            ODE units are not supported in steady-state EO mode).
    """

    def __init__(self, config: dict, backend: str | None = None) -> None:
        logger.info("Initializing EOFlowsheet")
        self.config = config
        # Resolve backend: explicit arg > simulation.backend > simulation.default_backend > "scipy"
        sim_cfg = config.get("simulation", {})
        self.backend = (
            backend
            or sim_cfg.get("backend")
            or sim_cfg.get("default_backend", "scipy")
        )
        self._unit_objects: dict = {}
        self._provider_map: dict = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self) -> dict:
        """
        Build and solve the EO system.

        After this method returns, ``self.x0`` holds the initial-guess vector
        and ``self.var_names`` holds a human-readable label for each element —
        both can be passed directly to
        :func:`processforge.provenance.build_run_info` for reproducibility.

        Returns:
            Stream result dict in the same format as ``Flowsheet.run()``:
            ``{stream_name: {"T": ..., "P": ..., "flowrate": ..., "z": {...}}}``
        """
        from processforge.providers.manager import teardown_providers
        manager = self._build()
        try:
            x0 = self._warm_start(manager)
            self.x0: "np.ndarray" = x0.copy()
            self.var_names: list[str] = self._build_var_names(manager)
            solver = EOSolver(backend=self.backend)
            x_sol, converged, stats = solver.solve(manager, x0)
            if not converged:
                logger.warning("EOFlowsheet: solver did not fully converge.")
            results = manager.extract_results(x_sol)
            logger.info("EOFlowsheet: simulation complete.")
            return results
        finally:
            teardown_providers(self._provider_map)

    # ------------------------------------------------------------------
    # Build phase
    # ------------------------------------------------------------------

    def _build(self) -> GlobalJacobianManager:
        """Parse config, build units, register everything with the manager."""
        self._check_no_tank()
        self._check_provider_backend_compat()

        from processforge.providers.manager import build_provider_map
        self._provider_map = build_provider_map(
            self.config.get("providers", {}), self.config
        )

        components = self._collect_components()
        logger.info(f"EOFlowsheet: components = {components}")

        manager = GlobalJacobianManager(components)

        # Register feed streams (boundary conditions)
        for name, stream_dict in self.config["streams"].items():
            manager.register_feed(name, stream_dict)

        # Build unit objects and register them
        self._unit_objects = self._build_units(self._provider_map)
        for unit_name, unit in self._unit_objects.items():
            unit_cfg = self.config["units"][unit_name]
            inlet_names = _get_inlets(unit_cfg)
            outlet_names = _get_outlets(unit_cfg)
            manager.register_unit(unit, inlet_names, outlet_names, unit_cfg)
            logger.debug(
                f"Registered unit '{unit_name}': "
                f"inlets={inlet_names} outlets={outlet_names}"
            )

        return manager

    def _build_units(self, provider_map: dict | None = None) -> dict:
        """Instantiate unit objects using the same factory logic as Flowsheet."""
        from processforge.units.pump import Pump
        from processforge.units.valve import Valve
        from processforge.units.strainer import Strainer
        from processforge.units.pipes import Pipes
        from processforge.units.heater import Heater
        from processforge.units.flash import Flash

        unit_types = {
            "Pump": (Pump, "kwargs"),
            "Valve": (Valve, "kwargs"),
            "Strainer": (Strainer, "kwargs"),
            "Pipes": (Pipes, "kwargs"),
            "Heater": (Heater, "params"),
            "Flash": (Flash, "params"),
        }

        from processforge.providers.manager import get_unit_provider
        _pm = provider_map or {}

        units = {}
        for unit_name, unit_cfg in self.config["units"].items():
            unit_type = unit_cfg["type"]
            if unit_type not in unit_types:
                raise ValueError(f"Unknown unit type '{unit_type}'.")
            cls, style = unit_types[unit_type]

            exclude = {"type", "in", "out", "provider"}
            if unit_type == "Flash":
                # Flash uses out_vap / out_liq — exclude both from kwargs
                exclude |= {"out_vap", "out_liq"}

            params = {k: v for k, v in unit_cfg.items() if k not in exclude}
            if style == "params":
                # Heater / Flash use __init__(name, params_dict)
                unit = cls(unit_name, params)
            else:
                unit = cls(unit_name, **params)

            # Attach provider so EO residuals can route thermo calls.
            if _pm:
                unit._provider = get_unit_provider(_pm, unit_cfg)

            units[unit_name] = unit
            logger.info(f"Built unit '{unit_name}' ({unit_type})")

        return units

    def _check_provider_backend_compat(self) -> None:
        """Raise if a Cantera provider is used with a symbolic EO backend.

        Cantera thermo calls return Python floats and are not symbolically
        differentiable.  They cannot be used with Pyomo or CasADi backends.
        """
        providers_cfg = self.config.get("providers", {})
        has_cantera = any(
            p.get("type") == "cantera" for p in providers_cfg.values()
        )
        if has_cantera and self.backend in ("pyomo", "casadi"):
            raise NotImplementedError(
                f"CanteraProvider is not compatible with backend='{self.backend}'. "
                "Cantera thermo calls are not symbolically differentiable. "
                "Use backend='scipy' (or omit 'backend' — it defaults to 'scipy')."
            )

    def _check_no_tank(self) -> None:
        """Raise if any Tank unit is present (dynamic units not EO-compatible)."""
        tank_units = [
            n for n, cfg in self.config["units"].items()
            if cfg.get("type") == "Tank"
        ]
        if tank_units:
            raise NotImplementedError(
                f"EOFlowsheet does not support Tank units: {tank_units}. "
                "Tank requires ODE/DAE integration. "
                "Use mode='dynamic' for flowsheets containing Tank units."
            )

    def _collect_components(self) -> list[str]:
        """
        Return sorted list of all component names found anywhere in the config.
        """
        comps: set[str] = set()
        for stream in self.config["streams"].values():
            comps.update(stream.get("z", {}).keys())
        return sorted(comps)

    # ------------------------------------------------------------------
    # Provenance helpers
    # ------------------------------------------------------------------

    def _build_var_names(self, manager: "GlobalJacobianManager") -> list[str]:
        """Return a human-readable label for every scalar in the x vector.

        Variable order mirrors ``GlobalJacobianManager`` layout:
        ``T``, ``P``, ``flowrate``, then one entry per component ``z_<comp>``
        for each registered stream in registration order.
        """
        suffixes = ["T", "P", "flowrate"] + [
            f"z_{c}" for c in manager.components
        ]
        names: list[str] = []
        for stream_name in manager._streams:
            for suffix in suffixes:
                names.append(f"{stream_name}/{suffix}")
        return names

    # ------------------------------------------------------------------
    # Warm-start
    # ------------------------------------------------------------------

    def _warm_start(self, manager: GlobalJacobianManager) -> "import numpy; numpy.ndarray":
        """
        Build initial guess x0 by one forward SM pass (no iteration).

        Feed stream values are propagated through the processing order
        (skipping back-edges / cycle streams).  Unresolved streams default
        to the first feed stream's values.
        """
        import numpy as np

        stream_vals: dict[str, dict] = {}

        # Seed with feed streams
        for name, feed in self.config["streams"].items():
            stream_vals[name] = deepcopy(feed)

        default_feed = deepcopy(next(iter(self.config["streams"].values())))

        # One forward pass in config order (no topological sort needed for warm start)
        for unit_name, unit_cfg in self.config["units"].items():
            unit = self._unit_objects[unit_name]
            inlet_names = _get_inlets(unit_cfg)
            outlet_names = _get_outlets(unit_cfg)

            # Skip if any inlet is not yet available
            if not all(n in stream_vals for n in inlet_names):
                continue

            inlet_snapshot = self._merge_inlets(stream_vals, inlet_names)

            try:
                if unit_cfg.get("type") == "Heater":
                    # Heater.run() takes full streams dict
                    tmp = {unit.inlet: inlet_snapshot}
                    result = unit.run(tmp)
                    stream_vals[outlet_names[0]] = result[unit.outlet]
                elif unit_cfg.get("type") == "Flash":
                    # Flash.run() takes full streams dict, produces 2 outlets
                    tmp = {unit.inlet: inlet_snapshot}
                    result = unit.run(tmp)
                    if len(outlet_names) >= 2:
                        stream_vals[outlet_names[0]] = result.get(unit.out_liq, {})
                        stream_vals[outlet_names[1]] = result.get(unit.out_vap, {})
                else:
                    out = unit.run(inlet_snapshot)
                    if outlet_names:
                        stream_vals[outlet_names[0]] = out
            except Exception as exc:  # noqa: BLE001
                logger.debug(
                    f"Warm-start: unit '{unit_name}' failed ({exc}); "
                    "using feed defaults."
                )

        # Fill any unresolved streams with defaults
        for name in manager._streams:
            if name not in stream_vals:
                stream_vals[name] = deepcopy(default_feed)

        return manager.get_x0(stream_vals)

    @staticmethod
    def _merge_inlets(stream_vals: dict, inlet_names: list[str]) -> dict:
        """Flow-weighted merge of multiple inlet stream dicts."""
        if len(inlet_names) == 1:
            return deepcopy(stream_vals[inlet_names[0]])

        merged: dict = {"T": 0.0, "P": 0.0, "flowrate": 0.0, "z": {}}
        total_f = 0.0
        for name in inlet_names:
            s = stream_vals.get(name, {})
            f = s.get("flowrate", 0.0)
            total_f += f
            merged["T"] += s.get("T", 298.15) * f
            merged["P"] += s.get("P", 101325.0) * f
            for c, frac in s.get("z", {}).items():
                merged["z"][c] = merged["z"].get(c, 0.0) + frac * f

        if total_f > 1e-20:
            merged["T"] /= total_f
            merged["P"] /= total_f
            for c in merged["z"]:
                merged["z"][c] /= total_f
        merged["flowrate"] = total_f
        return merged
