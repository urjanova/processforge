"""Native hydraulic network solver for pressure-driven flowsheets.

Builds a small nonlinear system from the flowsheet topology and solves it
with :func:`scipy.optimize.least_squares`:

* **Nodes** are streams.  Streams attached to a ``Tank`` or declared as feed
  streams are fixed-pressure boundaries; every other stream is a free node
  with a mass balance.
* **Edges** are single-inlet/single-outlet units — ``Pump``, ``Valve``,
  ``Pipes``, ``Strainer`` and (as a zero-resistance pass-through) ``Heater``.
  Each edge contributes one pressure–flow relation.
* **Tanks** are pressure boundaries with storage.  Their inlet and outlet
  flows are independent; the inventory is integrated by the flowsheet driver.

Check valves and pumps are non-return elements.  They are handled with an
active-set iteration: an element that wants reverse flow is closed (flow
fixed to zero, pressure relation dropped) and the network is re-solved.

The element pressure–flow laws live in
:mod:`processforge.units.hydraulics` so they stay independent of the
network assembly and can be unit-tested in isolation.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from loguru import logger
from scipy.optimize import least_squares

from .units.hydraulics import (
    HydraulicConfigError,
    HydraulicSolveError,
    design_point_delta_p,
    fit_pump_curve,
    pipe_delta_p,
    pump_curve_points_from_config,
    pump_pressure_rise,
    resolve_hydraulic_props,
    valve_delta_p,
)
from .utils.topology import get_inlets, get_outlets

# Unit types supported as hydraulic elements (single inlet, single outlet).
_HYDRAULIC_EDGE_TYPES = frozenset({"Pump", "Valve", "Pipes", "Strainer", "Heater"})
_PASS_THROUGH_TYPES = frozenset({"Heater"})


@dataclass
class HydraulicElement:
    """A single pressure–flow element (one unit, one inlet, one outlet)."""

    name: str
    kind: str  # "pump" | "valve" | "pipes" | "strainer" | "heater"
    inlet: str
    outlet: str
    params: dict
    nonreturn: bool = False


@dataclass
class HydraulicSolution:
    """Result of one hydraulic network solve."""

    flows: dict[str, float] = field(default_factory=dict)  # unit name -> F [mol/s]
    stream_flows: dict[str, float] = field(default_factory=dict)  # stream -> F [mol/s]
    pressures: dict[str, float] = field(default_factory=dict)  # stream -> P [Pa]
    closed: list[str] = field(default_factory=list)  # closed non-return units
    residual_norm: float = float("inf")
    iterations: int = 0

    def flow(self, unit_name: str) -> float:
        return self.flows.get(unit_name, 0.0)

    def pressure(self, stream_name: str) -> Optional[float]:
        return self.pressures.get(stream_name)


class HydraulicNetwork:
    """Pressure–flow network assembled from a processforge flowsheet."""

    def __init__(
        self,
        elements: list[HydraulicElement],
        fixed_pressures: dict[str, float],
        free_nodes: list[str],
        tol: float = 1e-8,
        max_iter: int = 50,
        flow_scale: float = 1e-3,
        pressure_scale: float = 1e5,
        node_aliases: dict[str, str] | None = None,
    ) -> None:
        self.elements = elements
        self.node_aliases = dict(node_aliases or {})
        self.tol = float(tol)
        self.max_iter = int(max_iter)
        self.flow_scale = float(flow_scale)
        self.pressure_scale = float(pressure_scale)
        self._last_x: Optional[np.ndarray] = None
        self._last_flows: dict[str, float] = {}

        # Every stream name that appears anywhere.  Used to expand the solved
        # canonical-node pressures back onto the individual streams.
        stream_names: set[str] = set(self.node_aliases) | set(self.node_aliases.values())
        for el in elements:
            stream_names.add(el.inlet)
            stream_names.add(el.outlet)
        stream_names.update(fixed_pressures)
        stream_names.update(free_nodes)
        self._stream_names = sorted(stream_names)

        # Merge aliased streams into a single canonical node and collapse the
        # fixed pressures onto those nodes.  Junction units alias their extra
        # inlets to the first inlet, so all of them share one node/unknown.
        canonical_fixed: dict[str, float] = {}
        for node, pressure in fixed_pressures.items():
            canon = self._canonical(node)
            if canon in canonical_fixed and not math.isclose(
                canonical_fixed[canon], float(pressure), rel_tol=1e-9, abs_tol=1e-6
            ):
                raise HydraulicConfigError(
                    "Junction stream pressure conflict: streams merged at a "
                    f"common node carry {canonical_fixed[canon]!r} Pa and "
                    f"{pressure!r} Pa."
                )
            canonical_fixed[canon] = float(pressure)
        self.fixed_pressures = canonical_fixed
        self.free_nodes = sorted(
            {self._canonical(n) for n in free_nodes} - set(canonical_fixed)
        )

        self._node_in: dict[str, list[int]] = {n: [] for n in self.free_nodes}
        self._node_out: dict[str, list[int]] = {n: [] for n in self.free_nodes}
        for idx, el in enumerate(self.elements):
            outlet = self._canonical(el.outlet)
            inlet = self._canonical(el.inlet)
            if outlet in self._node_in:
                self._node_in[outlet].append(idx)
            if inlet in self._node_out:
                self._node_out[inlet].append(idx)

        self._check_components_have_boundary()

    def _canonical(self, node: str) -> str:
        """Resolve a stream name to its junction-canonical node."""
        seen: set[str] = set()
        while node in self.node_aliases and node not in seen:
            seen.add(node)
            node = self.node_aliases[node]
        return node

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def solve(self) -> HydraulicSolution:
        """Solve the network, closing non-return elements as needed."""
        closed = [False] * len(self.elements)
        for i, el in enumerate(self.elements):
            if el.kind == "valve" and el.params.get("opening", 1.0) <= 0.0:
                closed[i] = True

        max_active = 10
        solution: Optional[HydraulicSolution] = None
        for _ in range(max_active):
            solution = self._solve_once(closed)
            newly_closed: list[int] = []
            for i, el in enumerate(self.elements):
                if not el.nonreturn or closed[i]:
                    continue
                f = solution.flows[el.name]
                p_in = solution.pressures[el.inlet]
                p_out = solution.pressures[el.outlet]
                flow_tol = self.tol * self.flow_scale
                pressure_tol = self.tol * self.pressure_scale
                if f < -flow_tol or (el.kind == "valve" and p_out > p_in + pressure_tol):
                    newly_closed.append(i)
            if not newly_closed:
                solution.closed = [
                    el.name for i, el in enumerate(self.elements) if closed[i]
                ]
                return solution
            for i in newly_closed:
                closed[i] = True
        raise HydraulicSolveError(
            "Hydraulic active-set iteration did not stabilize after "
            f"{max_active} rounds (closed elements: "
            f"{[el.name for i, el in enumerate(self.elements) if closed[i]]})"
        )

    def inlet_flow(self, stream_names: list[str]) -> float:
        """Total flow [mol/s] into a set of streams (used for tank inlets)."""
        targets = {self._canonical(s) for s in stream_names}
        total = 0.0
        for idx, el in enumerate(self.elements):
            if self._canonical(el.outlet) in targets:
                total += self._last_flows.get(el.name, 0.0)
        return total

    def outlet_flow(self, stream_names: list[str]) -> float:
        """Total flow [mol/s] out of a set of streams (used for tank outlets)."""
        targets = {self._canonical(s) for s in stream_names}
        total = 0.0
        seen: set[int] = set()
        for idx, el in enumerate(self.elements):
            if idx in seen:
                continue
            if self._canonical(el.inlet) in targets:
                seen.add(idx)
                total += self._last_flows.get(el.name, 0.0)
        return total

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    def _check_components_have_boundary(self) -> None:
        """Every connected component of stream nodes needs a pressure boundary."""
        adjacency: dict[str, set[str]] = {}
        for el in self.elements:
            inlet = self._canonical(el.inlet)
            outlet = self._canonical(el.outlet)
            adjacency.setdefault(inlet, set()).add(outlet)
            adjacency.setdefault(outlet, set()).add(inlet)
        for node in self.free_nodes:
            adjacency.setdefault(node, set())

        visited: set[str] = set()
        for start in adjacency:
            if start in visited:
                continue
            stack = [start]
            component: set[str] = set()
            has_boundary = start in self.fixed_pressures
            while stack:
                node = stack.pop()
                if node in component:
                    continue
                component.add(node)
                if node in self.fixed_pressures:
                    has_boundary = True
                for neighbour in adjacency.get(node, ()):
                    if neighbour not in component:
                        stack.append(neighbour)
            visited |= component
            if not has_boundary:
                raise HydraulicConfigError(
                    "Hydraulic network component "
                    f"{sorted(component)} has no pressure boundary. "
                    "Add a Tank or a feed stream to anchor the pressure."
                )

    def _initial_guess(self, closed: list[bool]) -> np.ndarray:
        n_free = len(self.free_nodes)
        n_edges = len(self.elements)
        x0 = np.zeros(n_free + n_edges)

        fixed_values = list(self.fixed_pressures.values())
        default_pressure = (
            float(np.mean(fixed_values)) if fixed_values else 101325.0
        )
        for i, node in enumerate(self.free_nodes):
            if self._last_x is not None and i < len(self._last_x):
                x0[i] = self._last_x[i]
            else:
                x0[i] = default_pressure

        # Pressure lookup including both fixed and (current guess) free nodes.
        pressures: dict[str, float] = {}
        for i, node in enumerate(self.free_nodes):
            pressures[node] = float(x0[i])
        for node, p in self.fixed_pressures.items():
            pressures[node] = float(p)

        for j, el in enumerate(self.elements):
            if self._last_x is not None and n_free + j < len(self._last_x):
                x0[n_free + j] = self._last_x[n_free + j]
            else:
                design = el.params.get("design_flowrate")
                if design and design > 0:
                    x0[n_free + j] = float(design)
                    continue
                # Estimate flow from the available pressure drop and the
                # element's resistance law so the solver starts near the
                # solution instead of at a generic 1e-3.
                p_in = pressures.get(self._canonical(el.inlet), default_pressure)
                p_out = pressures.get(self._canonical(el.outlet), default_pressure)
                dp = abs(p_in - p_out)
                if el.kind == "valve":
                    cv = el.params.get("cv", 1e-6)
                    f = el.params.get("opening", 1.0)
                    x0[n_free + j] = max(cv * f * math.sqrt(dp) if dp > 0 else 1e-4, 1e-6)
                elif el.kind == "pump":
                    # Start at a small positive flow; the pump curve determines
                    # the operating point.
                    x0[n_free + j] = 1e-4
                elif el.kind in ("pipes", "strainer"):
                    delta_p = el.params.get("delta_p", 1000.0)
                    design = el.params.get("design_flowrate", 1e-4)
                    x0[n_free + j] = design
                else:
                    x0[n_free + j] = 1e-4
        if closed:
            for j, is_closed in enumerate(closed):
                if is_closed:
                    x0[n_free + j] = 0.0
        return x0

    def _solve_once(self, closed: list[bool]) -> HydraulicSolution:
        n_free = len(self.free_nodes)
        n_edges = len(self.elements)

        # Short-circuit: nothing to solve (degenerate network).
        if n_free == 0 and n_edges == 0:
            return HydraulicSolution(
                flows={},
                stream_flows={},
                pressures={**self.fixed_pressures},
                residual_norm=0.0,
                iterations=0,
            )

        x0 = self._initial_guess(closed)

        x_scale = np.empty(n_free + n_edges)
        x_scale[:n_free] = self.pressure_scale
        x_scale[n_free:] = self.flow_scale

        result = least_squares(
            lambda x: self._residuals(x, closed),
            x0,
            x_scale=x_scale,
            ftol=self.tol,
            xtol=0.0,
            gtol=0.0,
            max_nfev=max(500, self.max_iter * 20),
        )
        self._last_x = result.x.copy()

        residual_norm = float(np.max(np.abs(result.fun))) if result.fun.size else 0.0
        if residual_norm > self.tol:
            raise HydraulicSolveError(
                self._failure_message(result, residual_norm, closed)
            )

        canonical_pressures: dict[str, float] = {}
        for i, node in enumerate(self.free_nodes):
            canonical_pressures[node] = float(result.x[i])
        canonical_pressures.update(self.fixed_pressures)

        # Expand canonical-node pressures onto every individual stream so a
        # junction's merged inlets all report the same pressure.
        pressures: dict[str, float] = {}
        for name in self._stream_names:
            pressures[name] = float(canonical_pressures[self._canonical(name)])

        flows: dict[str, float] = {}
        for j, el in enumerate(self.elements):
            flows[el.name] = float(result.x[n_free + j])
        self._last_flows = flows

        produced_by: dict[str, list[int]] = {}
        for j, el in enumerate(self.elements):
            produced_by.setdefault(el.outlet, []).append(j)
        stream_flows: dict[str, float] = {}
        for node in self._stream_names:
            if node in produced_by:
                stream_flows[node] = sum(
                    flows[self.elements[j].name] for j in produced_by[node]
                )
            else:
                # A declared feed (or merged boundary) node: report the total
                # flow leaving the node through its canonical edges.
                target = self._canonical(node)
                total = 0.0
                seen: set[int] = set()
                for j, el in enumerate(self.elements):
                    if j in seen:
                        continue
                    if self._canonical(el.inlet) == target:
                        seen.add(j)
                        total += flows[el.name]
                stream_flows[node] = total

        return HydraulicSolution(
            flows=flows,
            stream_flows=stream_flows,
            pressures=pressures,
            residual_norm=residual_norm,
            iterations=int(result.nfev),
        )

    def _equation_labels(self) -> list[str]:
        """Human-readable labels aligned with the residual vector."""
        labels = [f"mass balance at '{node}'" for node in self.free_nodes]
        labels.extend(
            f"{el.kind} '{el.name}' pressure-flow relation"
            for el in self.elements
        )
        return labels

    def _failure_message(
        self,
        result,
        residual_norm: float,
        closed: list[bool],
    ) -> str:
        """Build a diagnostic message for a non-converged hydraulic solve."""
        labels = self._equation_labels()
        worst = "unknown equation"
        if result.fun.size:
            worst_idx = int(np.argmax(np.abs(result.fun)))
            if 0 <= worst_idx < len(labels):
                worst = labels[worst_idx]
        closed_names = [
            el.name for i, el in enumerate(self.elements) if closed[i]
        ]
        active = (
            f"closed: {closed_names}" if closed_names else "all elements open"
        )
        return (
            "Hydraulic network did not converge "
            f"(scaled residual {residual_norm:.3e} > tol {self.tol:.3e}; "
            f"solver status {result.status}: {result.message}). "
            f"Worst equation: {worst}. Valve/pump active set: {active}. "
            "The flowsheet may be over-determined, missing a resistance element "
            "between pressure boundaries, or a boundary stream may be unanchored."
        )

    def _residuals(self, x: np.ndarray, closed: list[bool]) -> np.ndarray:
        n_free = len(self.free_nodes)
        pressures: dict[str, float] = {}
        for i, node in enumerate(self.free_nodes):
            pressures[node] = float(x[i])
        for node, p in self.fixed_pressures.items():
            pressures[node] = float(p)

        flows = x[n_free:]
        residuals = np.empty(n_free + len(self.elements))

        # --- node mass balances -----------------------------------------
        for i, node in enumerate(self.free_nodes):
            net = 0.0
            for j in self._node_in[node]:
                net += flows[j]
            for j in self._node_out[node]:
                net -= flows[j]
            residuals[i] = net / self.flow_scale

        # --- element pressure relations ---------------------------------
        for j, el in enumerate(self.elements):
            idx = n_free + j
            if closed[j]:
                residuals[idx] = flows[j] / self.flow_scale
                continue

            p_in = pressures[self._canonical(el.inlet)]
            p_out = pressures[self._canonical(el.outlet)]
            f = float(flows[j])

            if el.kind == "valve":
                dp = valve_delta_p(
                    f,
                    el.params["cv"],
                    el.params.get("opening", 1.0),
                    el.params.get("valve_function", "linear"),
                    el.params.get("alpha", 100.0),
                )
                residuals[idx] = (p_out - p_in + dp) / self.pressure_scale
            elif el.kind == "pipes":
                if el.params.get("mode") == "geometry":
                    dp = pipe_delta_p(
                        f,
                        el.params["diameter"],
                        el.params["length"],
                        el.params["viscosity"],
                        el.params["molar_volume"],
                        el.params["density"],
                        el.params.get("roughness", 0.0),
                    )
                else:
                    dp = design_point_delta_p(
                        f,
                        el.params["delta_p"],
                        el.params["design_flowrate"],
                    )
                residuals[idx] = (p_out - p_in + dp) / self.pressure_scale
            elif el.kind == "strainer":
                dp = design_point_delta_p(
                    f,
                    el.params["delta_p"],
                    el.params["design_flowrate"],
                )
                residuals[idx] = (p_out - p_in + dp) / self.pressure_scale
            elif el.kind == "pump":
                rise = pump_pressure_rise(
                    f,
                    el.params.get("delta_p"),
                    el.params.get("curve"),
                )
                residuals[idx] = (p_out - p_in - rise) / self.pressure_scale
            elif el.kind == "heater":
                residuals[idx] = (p_out - p_in) / self.pressure_scale
            else:
                raise HydraulicConfigError(
                    f"Unsupported hydraulic element kind '{el.kind}' on unit '{el.name}'"
                )
        return residuals


# ---------------------------------------------------------------------------
# Construction from a flowsheet config
# ---------------------------------------------------------------------------

def _preprocess_valve(cfg: dict) -> dict:
    cv = cfg.get("Cv", cfg.get("cv"))
    if cv is None or cv <= 0.0:
        raise HydraulicConfigError(
            "Valves in hydraulic mode need a positive IDAES-style 'Cv' "
            "[mol/s/Pa**0.5]. 'pressure_ratio' is not a flow model."
        )
    return {
        "cv": float(cv),
        "opening": float(cfg.get("valve_opening", cfg.get("opening", 1.0))),
        "valve_function": cfg.get("valve_function", "linear"),
        "alpha": float(cfg.get("alpha", 100.0)),
        "check_valve": bool(cfg.get("check_valve", False)),
    }


def _preprocess_pipes(cfg: dict, props: dict) -> dict:
    diameter = cfg.get("diameter")
    length = cfg.get("length", 1.0)
    viscosity = props.get("viscosity")
    molar_volume = props.get("molar_volume")
    density = props.get("density")
    geometry_ok = all(
        v is not None and v > 0.0
        for v in (diameter, length, viscosity, molar_volume, density)
    )
    if geometry_ok:
        return {
            "mode": "geometry",
            "diameter": float(diameter),
            "length": float(length),
            "viscosity": float(viscosity),
            "molar_volume": float(molar_volume),
            "density": float(density),
            "roughness": float(cfg.get("roughness", 0.0)),
        }
    delta_p = cfg.get("delta_p")
    design_flowrate = cfg.get("design_flowrate")
    if delta_p is not None and design_flowrate is not None and design_flowrate > 0.0:
        return {
            "mode": "design",
            "delta_p": float(delta_p),
            "design_flowrate": float(design_flowrate),
        }
    missing = []
    if viscosity is None:
        missing.append("viscosity")
    if molar_volume is None:
        missing.append("molar_volume/density")
    if diameter is None:
        missing.append("diameter")
    raise HydraulicConfigError(
        "Pipes in hydraulic mode need either full geometry "
        "(diameter, length, viscosity, molar_volume/density) or a design "
        f"point (delta_p + design_flowrate). Missing: {sorted(set(missing))}"
    )


def _preprocess_strainer(cfg: dict) -> dict:
    delta_p = cfg.get("delta_p")
    design_flowrate = cfg.get("design_flowrate")
    if delta_p is None or design_flowrate is None or design_flowrate <= 0.0:
        raise HydraulicConfigError(
            "Strainers in hydraulic mode need 'delta_p' at "
            "'design_flowrate' so a resistance can be derived."
        )
    return {
        "delta_p": float(delta_p),
        "design_flowrate": float(design_flowrate),
    }


def _preprocess_pump(cfg: dict, props: dict) -> dict:
    params: dict = {}
    curve_cfg = cfg.get("pump_curve")
    if curve_cfg is not None:
        density = props.get("density")
        points = pump_curve_points_from_config(curve_cfg, density)
        params["curve"] = fit_pump_curve(points)
    else:
        delta_p = cfg.get("delta_p")
        if delta_p is None or delta_p < 0.0:
            raise HydraulicConfigError(
                "Pumps in hydraulic mode need either 'delta_p' or 'pump_curve'"
            )
        params["delta_p"] = float(delta_p)
    if cfg.get("efficiency") is not None:
        params["efficiency"] = float(cfg["efficiency"])
    return params


def _preprocess_heater(_cfg: dict) -> dict:
    return {}


def build_hydraulic_network(
    config: dict,
    units: dict,
    initial_states: dict[str, dict],
    hydraulics_cfg: dict | None = None,
) -> tuple[HydraulicNetwork, dict[str, list[str]]]:
    """Build a :class:`HydraulicNetwork` from a flowsheet configuration.

    Args:
        config: The flowsheet configuration dict.
        units: Instantiated unit objects (from ``Flowsheet.build_units``).
        initial_states: Approximate stream states used to resolve provider
            properties for each element (from a forward pass at t=0).
        hydraulics_cfg: The ``simulation.hydraulics`` block.

    Returns:
        ``(network, tank_outlet_streams)`` where ``tank_outlet_streams`` maps
        each Tank unit name to the streams it feeds downstream.
    """
    hydraulics_cfg = hydraulics_cfg or {}

    fixed_pressures: dict[str, float] = {}
    for name, stream in config.get("streams", {}).items():
        p = stream.get("P")
        if p is None:
            raise HydraulicConfigError(f"Feed stream '{name}' needs a pressure")
        fixed_pressures[name] = float(p)

    tank_outlet_streams: dict[str, list[str]] = {}
    for name, cfg in config.get("units", {}).items():
        if cfg.get("type") != "Tank":
            continue
        p = cfg.get("P")
        if p is None:
            raise HydraulicConfigError(f"Tank '{name}' needs a pressure")
        connected = get_inlets(cfg) + get_outlets(cfg)
        for stream in connected:
            if stream in fixed_pressures and fixed_pressures[stream] != p:
                raise HydraulicConfigError(
                    f"Stream '{stream}' is anchored to conflicting pressures: "
                    f"{fixed_pressures[stream]} Pa and {p} Pa"
                )
            fixed_pressures[stream] = float(p)
        tank_outlet_streams[name] = list(get_outlets(cfg))

    # Streams that feed a downstream unit (Tank inlets included).  Anything a
    # unit produces but nobody consumes is a free discharge (drain/product)
    # and must be anchored at a pressure.
    consumed: set[str] = set()
    for cfg in config.get("units", {}).values():
        consumed.update(get_inlets(cfg))

    elements: list[HydraulicElement] = []
    node_aliases: dict[str, str] = {}
    for name, cfg in config.get("units", {}).items():
        unit_type = cfg.get("type", "")
        if unit_type in ("Tank", "SolverUnit"):
            continue
        if unit_type == "Flash":
            raise HydraulicConfigError(
                f"Unit '{name}' (Flash) is not supported in hydraulic mode"
            )
        if unit_type not in _HYDRAULIC_EDGE_TYPES:
            raise HydraulicConfigError(
                f"Unit '{name}' of type '{unit_type}' is not supported in "
                "hydraulic mode. Supported types: "
                f"{sorted(_HYDRAULIC_EDGE_TYPES | {'Tank'})}"
            )

        inlets = get_inlets(cfg)
        outlets = get_outlets(cfg)
        if not inlets:
            raise HydraulicConfigError(
                f"Unit '{name}' needs at least one inlet in hydraulic mode"
            )
        if len(outlets) != 1:
            raise HydraulicConfigError(
                f"Unit '{name}' must have exactly one outlet in hydraulic mode"
            )
        if unit_type == "Pump" and len(inlets) != 1:
            raise HydraulicConfigError(
                f"Pump '{name}' must have exactly one inlet in hydraulic mode"
            )

        inlet, outlet = inlets[0], outlets[0]
        # Multiple inlets are a zero-volume junction: merge the extra inlet
        # streams onto the first so the junction has one pressure unknown and
        # the node mass balance sums every incoming branch.
        for extra_inlet in inlets[1:]:
            if extra_inlet != inlet:
                node_aliases[extra_inlet] = inlet
        if len(inlets) > 1:
            logger.info(
                f"Unit '{name}': merging inlets {inlets[1:]} into junction "
                f"stream '{inlet}'"
            )

        # A produced stream that nothing consumes is a free discharge.  Anchor
        # it at an explicit ``outlet_pressure``, a declared stream pressure, or
        # atmospheric pressure.
        if outlet not in consumed and outlet not in fixed_pressures:
            discharge_p = cfg.get("outlet_pressure")
            if discharge_p is None:
                declared = config.get("streams", {}).get(outlet) or {}
                discharge_p = declared.get("P")
            if discharge_p is None:
                discharge_p = 101325.0
                logger.warning(
                    f"Stream '{outlet}' from unit '{name}' is a free discharge "
                    "with no 'outlet_pressure'; assuming atmospheric pressure "
                    "(101325 Pa)."
                )
            fixed_pressures[outlet] = float(discharge_p)
            logger.info(
                f"Stream '{outlet}' from unit '{name}' is a terminal outlet; "
                f"anchored at {float(discharge_p):.0f} Pa"
            )

        state = initial_states.get(inlet) or {}
        props = resolve_hydraulic_props(units[name], state)

        if unit_type == "Valve":
            if cfg.get("pressure_ratio") is not None:
                logger.warning(
                    f"Valve '{name}': 'pressure_ratio' is legacy and ignored in "
                    "hydraulic mode; the flow is governed by 'Cv'."
                )
            params = _preprocess_valve(cfg)
            nonreturn = bool(params["check_valve"])
        elif unit_type == "Pipes":
            params = _preprocess_pipes(cfg, props)
            nonreturn = False
        elif unit_type == "Strainer":
            params = _preprocess_strainer(cfg)
            nonreturn = False
        elif unit_type == "Pump":
            params = _preprocess_pump(cfg, props)
            nonreturn = True
        else:  # Heater
            params = _preprocess_heater(cfg)
            nonreturn = False

        elements.append(
            HydraulicElement(
                name=name,
                kind=unit_type.lower(),
                inlet=inlet,
                outlet=outlet,
                params=params,
                nonreturn=nonreturn,
            )
        )

    def _canonical(node: str) -> str:
        seen: set[str] = set()
        while node in node_aliases and node not in seen:
            seen.add(node)
            node = node_aliases[node]
        return node

    all_nodes: set[str] = set(fixed_pressures)
    for el in elements:
        all_nodes.add(_canonical(el.inlet))
        all_nodes.add(_canonical(el.outlet))
    free_nodes = sorted(all_nodes - {_canonical(n) for n in fixed_pressures})

    # Pick a flow scale that matches the flowsheet instead of the small water
    # demo default.  An explicit ``flow_scale`` always wins.
    if "flow_scale" in hydraulics_cfg:
        flow_scale = float(hydraulics_cfg["flow_scale"])
    else:
        candidates = [
            float(s.get("flowrate", 0.0))
            for s in config.get("streams", {}).values()
            if s.get("flowrate")
        ]
        candidates.extend(
            float(cfg["design_flowrate"])
            for cfg in config.get("units", {}).values()
            if cfg.get("design_flowrate")
        )
        flow_scale = max(candidates, default=1e-3)

    network = HydraulicNetwork(
        elements=elements,
        fixed_pressures=fixed_pressures,
        free_nodes=free_nodes,
        tol=float(hydraulics_cfg.get("tol", 1e-8)),
        max_iter=int(hydraulics_cfg.get("max_iter", 50)),
        flow_scale=flow_scale,
        pressure_scale=float(hydraulics_cfg.get("pressure_scale", 1e5)),
        node_aliases=node_aliases,
    )
    logger.info(
        f"Hydraulic network built: {len(elements)} elements, "
        f"{len(network.free_nodes)} free nodes, "
        f"{len(network.fixed_pressures)} boundaries"
    )
    return network, tank_outlet_streams
