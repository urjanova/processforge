"""EOFlowsheet: equation-oriented steady-state flowsheet solver."""
from __future__ import annotations

from copy import deepcopy

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

    def __init__(self, config: dict, backend: str = "pyomo") -> None:
        logger.info("Initializing EOFlowsheet")
        self.config = config
        self.backend = backend
        self._unit_objects: dict = {}

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self) -> dict:
        """
        Build and solve the EO system.

        Returns:
            Stream result dict in the same format as ``Flowsheet.run()``:
            ``{stream_name: {"T": ..., "P": ..., "flowrate": ..., "z": {...}}}``
        """
        manager = self._build()
        x0 = self._warm_start(manager)
        solver = EOSolver(backend=self.backend)
        x_sol, converged, stats = solver.solve(manager, x0)
        if not converged:
            logger.warning("EOFlowsheet: solver did not fully converge.")
        results = manager.extract_results(x_sol)
        logger.info("EOFlowsheet: simulation complete.")
        return results

    # ------------------------------------------------------------------
    # Build phase
    # ------------------------------------------------------------------

    def _build(self) -> GlobalJacobianManager:
        """Parse config, build units, register everything with the manager."""
        self._check_no_tank()

        components = self._collect_components()
        logger.info(f"EOFlowsheet: components = {components}")

        manager = GlobalJacobianManager(components)

        # Register feed streams (boundary conditions)
        for name, stream_dict in self.config["streams"].items():
            manager.register_feed(name, stream_dict)

        # Build unit objects and register them
        self._unit_objects = self._build_units()
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

    def _build_units(self) -> dict:
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

        units = {}
        for unit_name, unit_cfg in self.config["units"].items():
            unit_type = unit_cfg["type"]
            if unit_type not in unit_types:
                raise ValueError(f"Unknown unit type '{unit_type}'.")
            cls, style = unit_types[unit_type]

            exclude = {"type", "in", "out"}
            if unit_type == "Flash":
                # Flash uses out_vap / out_liq — exclude both from kwargs
                exclude |= {"out_vap", "out_liq"}

            params = {k: v for k, v in unit_cfg.items() if k not in exclude}
            if style == "params":
                # Heater / Flash use __init__(name, params_dict)
                units[unit_name] = cls(unit_name, params)
            else:
                units[unit_name] = cls(unit_name, **params)

            logger.info(f"Built unit '{unit_name}' ({unit_type})")

        return units

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
