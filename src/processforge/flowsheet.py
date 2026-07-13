from __future__ import annotations

from typing import TYPE_CHECKING

from loguru import logger
import numpy as np
from copy import deepcopy

if TYPE_CHECKING:
    from .providers.manager import ProviderMap

from .providers.manager import build_provider_map, teardown_providers
from .types import FlowsheetConfig, MergedInletTimeseries, StreamTimeseries
from .solver import Solver
from .providers.manager import UnitProviderConfig
from .units.registry import get_unit_class

_UNIT_META_KEYS = {"type", "in", "out", "material", "material_mix", "provider"}
_DEFAULT_T = 298.15
_DEFAULT_P = 101325.0


def _build_unit(unit_name, unit_config, provider_map):
    """Instantiate a single unit from its config dict."""
    unit_type = unit_config["type"]
    unit_class = get_unit_class(unit_type)
    if not unit_class:
        logger.error(f"Unknown unit type: {unit_type}")
        raise ValueError(f"Unknown unit type: {unit_type}")
    unit = unit_class(
        unit_name,
        **{k: v for k, v in unit_config.items() if k not in _UNIT_META_KEYS},
    )
    if "material" in unit_config:
        unit.material = unit_config["material"]
    if "material_mix" in unit_config:
        unit.material_mix = unit_config["material_mix"]
    if provider_map:
        unit._provider = provider_map.resolve(UnitProviderConfig(**unit_config))
    return unit


def _merge_provider_types(provider_map):
    """Register any extra unit types contributed by active providers."""
    if not provider_map:
        return
    for provider in provider_map.values():
        extra = getattr(provider, "unit_types", {})
        if extra:
            from .units.registry import register_unit
            for name, cls in extra.items():
                register_unit(name, cls)


class Flowsheet:
    """
    A class to represent and simulate a process flowsheet.
    The Flowsheet class manages units and streams for process simulations,
    supporting both steady-state and dynamic modes. It builds units from
    a configuration, runs simulations, and stores results.
    Attributes:
        config (dict): Configuration dictionary containing simulation settings,
                       units, and streams.
        units (dict): Dictionary of instantiated unit objects, keyed by name.
        streams (dict): Dictionary of stream data (initially empty).
        results (dict): Dictionary of simulation results, populated after running.
    """

    """
    Initializes the Flowsheet with the given configuration.
    Args:
        config (dict): Configuration dictionary including units, streams,
                        and simulation settings.
    """
    """
    Builds unit objects from the configuration.
    Iterates through the units in the config, instantiates the appropriate
    unit classes (e.g., Pump, Valve, etc.), and stores them in self.units.
    Raises:
        ValueError: If an unknown unit type is encountered.
    """
    """
    Starts the simulation based on the mode specified in the config.
    Builds units, determines the simulation mode (steady or dynamic),
    and calls the appropriate run method.
    Returns:
        dict: Simulation results.
    """
    """
    Runs a steady-state simulation.
    Processes each unit sequentially, updating stream results based on
    inlet conditions.
    Returns:
        dict: Dictionary of stream results after steady-state simulation.
    """
    """
    Runs a dynamic simulation.
    Simulates over time for dynamic units like Tank, using a solver.
    Assumes a simple 1-tank system for now.
    Args:
        solver: The solver object used for dynamic simulation.
    Returns:
        dict: Dictionary of stream results after dynamic simulation,
                including time-series data for dynamic units.
    """

    def __init__(self, config):
        logger.info("Initializing Flowsheet")
        self.config = config
        self.units = {}
        self.streams = {}
        self.tear_streams = set()
        self.has_recycles = False

    def build_units(self, provider_map: ProviderMap | None = None):
        logger.info("Building units")
        _merge_provider_types(provider_map=provider_map)
        self.units = {
            name: _build_unit(name, cfg, provider_map)
            for name, cfg in self.config["units"].items()
        }
        for name, unit in self.units.items():
            logger.info(f"Built unit {name} of type {self.config['units'][name]['type']}")

    def _get_unit_inlets(self, unit_name):
        """Get inlet stream(s) for a unit, handling both single and multiple inlets."""
        inlet = self.config["units"][unit_name].get("in")
        if isinstance(inlet, list):
            return inlet
        return [inlet] if inlet else []

    def _get_unit_outlets(self, unit_name):
        """Get all output streams for a unit."""
        cfg = self.config["units"][unit_name]
        outlets = []
        if "out" in cfg:
            out = cfg["out"]
            if isinstance(out, list):
                outlets.extend(out)
            else:
                outlets.append(out)
        if "retentate_out" in cfg:
            outlets.append(cfg["retentate_out"])
        if "permeate_out" in cfg:
            outlets.append(cfg["permeate_out"])
        # Add other potential outlet keys here
        return outlets

    def _detect_cycles_and_tear_streams(self):
        """
        Detects cycles in the flowsheet graph and identifies tear streams.
        Uses DFS to find back-edges (cycles). For each cycle found, selects a tear stream,
        preferring streams that enter Tank units (since Tanks provide holdup for stability).
        
        Returns:
            tuple: (has_cycles: bool, tear_streams: set, cycles: list of lists)
        """
        # Build graph: stream -> producing unit, unit -> output stream
        stream_producers = {}
        for unit_name in self.config["units"]:
            for out_stream in self._get_unit_outlets(unit_name):
                stream_producers[out_stream] = unit_name
        
        # Build adjacency: unit -> list of downstream units
        adj = {unit_name: [] for unit_name in self.units}
        for unit_name, unit_config in self.config["units"].items():
            inlets = self._get_unit_inlets(unit_name)
            for inlet_stream in inlets:
                if inlet_stream in stream_producers:
                    producer = stream_producers[inlet_stream]
                    if producer != unit_name:  # Avoid self-loops
                        adj[producer].append((unit_name, inlet_stream))
        
        # DFS to find cycles
        WHITE, GRAY, BLACK = 0, 1, 2
        color = {u: WHITE for u in self.units}
        parent = {u: None for u in self.units}
        cycles = []
        back_edges = []  # (from_unit, to_unit, stream_name)
        
        def dfs(u, path):
            color[u] = GRAY
            path.append(u)
            for v, stream in adj[u]:
                if color[v] == GRAY:
                    # Found a back-edge (cycle)
                    cycle_start = path.index(v)
                    cycle = path[cycle_start:] + [v]
                    cycles.append(cycle)
                    back_edges.append((u, v, stream))
                elif color[v] == WHITE:
                    parent[v] = u
                    dfs(v, path)
            path.pop()
            color[u] = BLACK
        
        for unit in self.units:
            if color[unit] == WHITE:
                dfs(unit, [])
        
        if not cycles:
            return False, set(), []
        
        # Select tear streams - prefer streams entering Tanks
        tear_streams = set()
        for from_unit, to_unit, stream in back_edges:
            tear_streams.add(stream)
            logger.info(f"Identified tear stream: '{stream}' (from {from_unit} to {to_unit})")
        
        return True, tear_streams, cycles

    def run(self):
        """
        In steady-state mode, the simulation models the system at equilibrium, where process variables (such as flows, pressures, and temperatures) are assumed to be constant over time, focusing on the long-term behavior without considering transients or time-dependent changes.
        In dynamic mode, the simulation accounts for time-dependent changes, modeling how the system evolves over time, including transients, start-ups, shut-downs, and responses to disturbances, providing a more comprehensive view of process dynamics.
        Runs the simulation based on the configured mode. This method initializes the simulation by building the units, determines the simulation mode from the configuration (defaulting to 'steady' if not specified), and executes either a dynamic or steady-state simulation accordingly.
        Returns:
            The result of the simulation run, which depends on the mode:
            - For dynamic mode: the output of self.run_dynamic(solver).
            - For steady mode: the output of self.run_steady().

        Example of difference between modes:


        Consider a tank filling with liquid from a pump and discharging through an outlet. In steady-state mode, the simulation calculates the final equilibrium level and flow rates once the system stabilizes, ignoring the time it takes to fill and the transients. In dynamic mode, it simulates the level rising and falling over time, showing how the tank level changes gradually in response to inflow and outflow until it reaches a steady value, capturing transients like initial surges or responses to flow changes.
        """

        logger.info("Starting simulation")
        flowsheet_cfg = FlowsheetConfig.from_dict(self.config)
        provider_map = build_provider_map(providers_config=flowsheet_cfg.providers, flowsheet_config=flowsheet_cfg)
        try:
            self.build_units(provider_map=provider_map)
            mode = self.config.get("simulation", {}).get("mode", "steady")
            logger.info(f"Simulation mode: {mode}")
            solver = Solver()
            return self.run_dynamic(solver) if mode == "dynamic" else self.run_steady()
        finally:
            teardown_providers(provider_map)

    def run_steady(self):
        """
        Runs a steady-state simulation for the flowsheet.
        This method initializes the results dictionary with the configured streams,
        then iterates through each unit in the flowsheet. For flowsheets with recycles,
        it uses iterative convergence with optional Wegstein acceleration.
        Returns:
            dict: A dictionary containing the updated stream results after simulation.
        """
        logger.info("Running steady-state simulation")
        
        # Determine processing order (this also detects recycles)
        processing_order = self._get_processing_order()
        
        if self.has_recycles:
            return self._run_steady_with_recycle(processing_order)
        
        # Simple sequential processing for non-recycle flowsheets
        results = {k: deepcopy(v) for k, v in self.config["streams"].items()}
        for unit_name in processing_order:
            unit = self.units[unit_name]
            cfg = self.config["units"][unit_name]
            inlet = self._get_merged_inlet(results, cfg)
            unit_out = unit.run(inlet)
            if "out" in cfg:
                out_val = cfg["out"]
                if isinstance(out_val, list):
                    for stream_name in out_val:
                        results[stream_name] = unit_out
                else:
                    results[out_val] = unit_out
            else:
                for out_key in ["retentate_out", "permeate_out"]:
                    if out_key in cfg and out_key in unit_out:
                        results[cfg[out_key]] = unit_out[out_key]
            logger.debug(f"Processed unit {unit_name}")
        self.results = results
        logger.info("Steady-state simulation completed")
        return results

    @staticmethod
    def _flow_weighted_merge_streams(streams: list[dict]) -> dict:
        """Flow-weighted merge of single-timestep stream dicts into one.

        Computes total flowrate as the sum, and flow-weighted averages for
        temperature, pressure, and mole fractions.
        """
        merged: dict = {"T": 0.0, "P": 0.0, "flowrate": 0.0, "z": {}}
        total_flow = 0.0

        for stream in streams:
            flow = stream.get("flowrate", 0.0)
            total_flow += flow
            merged["T"] += stream.get("T", _DEFAULT_T) * flow
            merged["P"] += stream.get("P", _DEFAULT_P) * flow
            for comp, frac in stream.get("z", {}).items():
                merged["z"].setdefault(comp, 0.0)
                merged["z"][comp] += frac * flow

        if total_flow > 0:
            merged["T"] /= total_flow
            merged["P"] /= total_flow
            for comp in merged["z"]:
                merged["z"][comp] /= total_flow

        merged["flowrate"] = total_flow
        return merged

    def _get_merged_inlet(self, results: dict, unit_cfg: dict) -> dict:
        """Get the inlet stream for a unit, merging multiple inlets if needed."""
        inlets = unit_cfg.get("in")
        if not inlets:
            return {}
        if isinstance(inlets, str):
            return results.get(inlets, {})
        return self._flow_weighted_merge_streams(
            [results.get(name, {}) for name in inlets]
        )

    def _run_steady_with_recycle(self, processing_order):
        """
        Runs steady-state simulation with recycle loop convergence.
        Uses iterative sequential modular approach with Wegstein acceleration.
        
        Args:
            processing_order: List of unit names in processing order.
            
        Returns:
            dict: Converged stream results.
        """
        logger.info("Running steady-state simulation with recycle convergence")
        
        max_iter = 100
        tolerance = 1e-6

        # Initialize results with feed streams
        results = {k: deepcopy(v) for k, v in self.config["streams"].items()}

        # Auto-initialize tear streams from the first available feed stream
        default_stream = deepcopy(next(iter(self.config["streams"].values())))
        for stream_name in self.tear_streams:
            if stream_name not in results:
                results[stream_name] = deepcopy(default_stream)
                logger.debug(f"Auto-initialized tear stream '{stream_name}' from feed values")
        
        # Wegstein acceleration storage
        wegstein_prev = {}  # Previous iteration values for each tear stream property
        wegstein_prev2 = {}  # Two iterations ago
        
        converged = False
        for iteration in range(max_iter):
            # Store old tear stream values for convergence check.
            # Manual copy avoids deepcopy overhead on large composition dicts.
            old_tear_values = {}
            for stream in self.tear_streams:
                if stream in results:
                    s = results[stream]
                    old_tear_values[stream] = {
                        "T": s["T"],
                        "P": s["P"],
                        "flowrate": s["flowrate"],
                        "z": dict(s.get("z", {})),
                    }
            
            # Process all units in order
            for unit_name in processing_order:
                unit = self.units[unit_name]
                cfg = self.config["units"][unit_name]
                inlet = self._get_merged_inlet(results, cfg)
                unit_out = unit.run(inlet)
                if "out" in cfg:
                    out_val = cfg["out"]
                    if isinstance(out_val, list):
                        for stream_name in out_val:
                            results[stream_name] = unit_out
                    else:
                        results[out_val] = unit_out
                else:
                    for out_key in ["retentate_out", "permeate_out"]:
                        if out_key in cfg and out_key in unit_out:
                            results[cfg[out_key]] = unit_out[out_key]
            
            # Check convergence on tear streams
            max_error = 0.0
            for stream_name in self.tear_streams:
                if stream_name not in results or stream_name not in old_tear_values:
                    continue
                    
                new_val = results[stream_name]
                old_val = old_tear_values[stream_name]
                
                # Compare flowrate, T, P
                for prop in ["flowrate", "T", "P"]:
                    if prop in new_val and prop in old_val:
                        old_p = old_val[prop]
                        new_p = new_val[prop]
                        if abs(old_p) > 1e-10:
                            error = abs(new_p - old_p) / abs(old_p)
                        else:
                            error = abs(new_p - old_p)
                        max_error = max(max_error, error)
                
                # Compare compositions
                if "z" in new_val and "z" in old_val:
                    for comp in new_val["z"]:
                        if comp in old_val["z"]:
                            old_z = old_val["z"][comp]
                            new_z = new_val["z"][comp]
                            if abs(old_z) > 1e-10:
                                error = abs(new_z - old_z) / abs(old_z)
                            else:
                                error = abs(new_z - old_z)
                            max_error = max(max_error, error)
            
            logger.debug(f"Iteration {iteration + 1}: max relative error = {max_error:.2e}")
            
            if max_error < tolerance:
                converged = True
                logger.info(f"Recycle converged after {iteration + 1} iterations (error: {max_error:.2e})")
                break
            
            # Apply Wegstein acceleration if we have enough history
            if iteration >= 2:
                for stream_name in self.tear_streams:
                    if stream_name in results:
                        results[stream_name] = self._apply_wegstein(
                            results[stream_name],
                            wegstein_prev.get(stream_name, {}),
                            wegstein_prev2.get(stream_name, {})
                        )
            
            # Update history — Wegstein only reads T, P, flowrate, so skip z.
            wegstein_prev2 = {k: {"T": v["T"], "P": v["P"], "flowrate": v["flowrate"]} for k, v in wegstein_prev.items()}
            wegstein_prev = {k: {"T": v["T"], "P": v["P"], "flowrate": v["flowrate"]} for k, v in old_tear_values.items()}
        
        if not converged:
            logger.warning(f"Recycle did not converge after {max_iter} iterations (error: {max_error:.2e})")
        
        self.results = results
        logger.info("Steady-state simulation with recycle completed")
        return results

    def _apply_wegstein(self, new_val, prev_val, prev2_val):
        """
        Apply Wegstein acceleration to speed up convergence.
        
        Wegstein formula: x_new = x + q * (g(x) - x)
        where q = s / (s - 1), s = (g(x) - g(x_prev)) / (x - x_prev)
        
        Bounded q to [-5, 0] for stability.
        """
        if not prev_val or not prev2_val:
            return new_val
        
        accelerated = {
            "T": new_val["T"],
            "P": new_val["P"],
            "flowrate": new_val["flowrate"],
            "z": dict(new_val.get("z", {})),
        }
        
        for prop in ["flowrate", "T", "P"]:
            if prop not in new_val or prop not in prev_val or prop not in prev2_val:
                continue
            
            x = prev_val.get(prop, 0)
            x_prev = prev2_val.get(prop, 0)
            g_x = new_val.get(prop, 0)
            g_x_prev = prev_val.get(prop, 0)
            
            dx = x - x_prev
            dg = g_x - g_x_prev
            
            if abs(dx) > 1e-12 and abs(dg - dx) > 1e-12:
                s = dg / dx
                q = s / (s - 1) if abs(s - 1) > 1e-12 else 0
                q = max(-5, min(0, q))  # Bound for stability
                accelerated[prop] = x + q * (g_x - x)
        
        return accelerated

    def _get_processing_order(self):
        """
        Determines the processing order of units in the flowsheet using topological sorting.
        This method constructs a directed graph where units are nodes and edges represent
        the flow of streams from producer units to consumer units. It then applies Kahn's
        algorithm (a breadth-first topological sort) to compute a valid processing order.
        
        For flowsheets with recycles (cycles), it identifies tear streams and computes
        a valid processing order that starts from units fed by external feeds or tear streams.
        
        Returns:
            list: A list of unit names in the order they should be processed.
        Raises:
            ValueError: If recycle loops exist without Tank units.
        """
        # First, detect cycles and tear streams
        has_cycles, tear_streams, cycles = self._detect_cycles_and_tear_streams()
        self.has_recycles = has_cycles
        self.tear_streams = tear_streams
        
        if has_cycles:
            logger.info(f"Detected {len(cycles)} recycle loop(s) with tear streams: {tear_streams}")

        stream_producers = {}
        for unit_name in self.config["units"]:
            for out_stream in self._get_unit_outlets(unit_name):
                stream_producers[out_stream] = unit_name

        # Create adjacency list (unit -> downstream units)
        adj = {unit_name: [] for unit_name in self.units}
        in_degree = dict.fromkeys(self.units, 0)

        for unit_name, unit_config in self.config["units"].items():
            inlets = self._get_unit_inlets(unit_name)
            for inlet_stream in inlets:
                # Skip tear streams when computing in-degree (break the cycle)
                if inlet_stream in tear_streams:
                    continue
                if inlet_stream in stream_producers:
                    producer_unit = stream_producers[inlet_stream]
                    adj[producer_unit].append(unit_name)
                    in_degree[unit_name] += 1

        # Kahn's algorithm for topological sort
        queue = [u for u, d in in_degree.items() if d == 0]
        order = []
        while queue:
            u = queue.pop(0)
            order.append(u)
            for v in adj.get(u, []):
                in_degree[v] -= 1
                if in_degree[v] == 0:
                    queue.append(v)

        if len(order) != len(self.units):
            # This shouldn't happen if tear streams are correctly identified
            remaining = [u for u in self.units if u not in order]
            logger.error(
                f"Could not determine processing order. Remaining units: {remaining}"
            )
            raise ValueError("Could not determine processing order for flowsheet.")

        logger.info(f"Processing order: {order}")
        return order

    def _process_static_unit(self, unit, inlet_stream_ts, t_eval, num_steps):
        """
        Helper method to process a static unit in dynamic simulation.
        Builds timeseries for the outlet stream by applying the unit's logic at each time step.
        Returns:
            dict: Timeseries dictionary for the outlet stream.
        """
        outlet_stream_ts = {prop: [] for prop in inlet_stream_ts if prop != "z"}
        outlet_stream_ts["time"] = t_eval.tolist()
        if "z" in inlet_stream_ts:
            outlet_stream_ts["z"] = {comp: [] for comp in inlet_stream_ts["z"].keys()}

        for i in range(num_steps):
            inlet_snapshot = {}
            for prop, values in inlet_stream_ts.items():

                if prop == "time":
                    continue
                if prop == "z":
                    inlet_snapshot["z"] = {
                        comp: comp_ts[i] for comp, comp_ts in values.items()
                    }
                    logger.debug(
                        f"Inlet snapshot for {prop} at step {i}: {inlet_snapshot['z']}"
                    )
                else:
                    if isinstance(values, list) and i < len(values):
                        logger.debug(
                            f"Inlet snapshot for {prop} at step {i}: {values[i]}"
                        )
                        inlet_snapshot[prop] = values[i]
                    else:
                        # For non-list properties, pass them as is.
                        inlet_snapshot[prop] = values

            outlet_snapshot = unit.run(inlet_snapshot)

            for prop, value in outlet_snapshot.items():
                logger.debug(f"Outlet snapshot for {prop}: {value}")
                if prop == "z":
                    for comp, comp_val in value.items():
                        if comp not in outlet_stream_ts["z"]:
                            outlet_stream_ts["z"][comp] = [0.0] * i
                        outlet_stream_ts["z"][comp].append(comp_val)
                else:
                    if prop not in outlet_stream_ts:
                        outlet_stream_ts[prop] = []
                    outlet_stream_ts[prop].append(value)

        for prop, values in outlet_stream_ts.items():
            if prop == "time":
                continue
            if prop == "z":
                for comp, comp_values in values.items():
                    if len(comp_values) < num_steps:
                        inlet_comp_ts = inlet_stream_ts.get("z", {}).get(
                            comp, [0.0] * num_steps
                        )
                        outlet_stream_ts["z"][comp] = inlet_comp_ts
            elif len(values) < num_steps:
                outlet_stream_ts[prop] = inlet_stream_ts.get(prop, [0.0] * num_steps)

        return outlet_stream_ts

    def run_dynamic(self, solver):
        """
        Run a dynamic simulation of the flowsheet over a specified time range.
        This method performs a time-dependent simulation by processing units in topological order.
        For units with a specialized 'run_dynamic' method, it delegates to that method for full
        timeseries handling. For static units, it applies their logic at each discrete time step,
        building timeseries data for inlet and outlet streams.
        
        For flowsheets with recycles, it uses timestep-by-timestep processing where recycle
        stream values from the previous timestep are used as inputs for the current timestep.
        The Tank holdup naturally integrates and stabilizes the recycle.
        
        Args:
            solver: The solver object (e.g., from scipy.integrate) used for dynamic unit simulations.
                    This is passed to unit-specific 'run_dynamic' methods if available.
        Returns:
            dict: A dictionary containing timeseries results for all streams in the flowsheet.
                  Each stream's data includes a 'time' key with the evaluation times, and other
                  properties as lists or nested dictionaries (e.g., 'z' for compositions).
        """
        logger.info("Running dynamic simulation")
        sim_config = self.config.get("simulation", {})
        start_time, end_time, time_step = (
            sim_config.get("t0", 0),
            sim_config.get("tf", 100),
            sim_config.get("dt", 1),
        )
        # Generate time evaluation points from start to end with given step size
        t_eval = np.arange(start_time, end_time + time_step, time_step)
        num_steps = len(t_eval)

        # Get processing order (this also detects recycles)
        processing_order = self._get_processing_order()

        # Initialize results with timeseries structure for all feed streams
        results = {}
        for stream_name, stream_data in self.config["streams"].items():
            results[stream_name] = self._init_stream_timeseries(stream_name, stream_data, t_eval, num_steps)

        if self.has_recycles:
            # Timestep-by-timestep processing for recycle loops
            return self._run_dynamic_with_recycle(solver, processing_order, results, t_eval, num_steps, start_time, end_time)
        else:
            # Standard unit-by-unit processing for non-recycle flowsheets
            return self._run_dynamic_sequential(solver, processing_order, results, t_eval, num_steps, start_time, end_time)

    def _init_stream_timeseries(self, stream_name, stream_data, t_eval, num_steps):
        """Initialize a stream's timeseries structure from its configuration."""
        stream_ts = {"time": t_eval.tolist()}
        for k, v in stream_data.items():
            if k == "z":
                if isinstance(v, dict):
                    stream_ts[k] = {
                        comp: [float(comp_val)] * num_steps
                        for comp, comp_val in v.items()
                    }
                else:
                    raise ValueError(
                        f"Stream '{stream_name}' property 'z' must be a dict of compositions."
                    )
            else:
                try:
                    stream_ts[k] = [float(v)] * num_steps
                except (ValueError, TypeError) as exc:
                    raise ValueError(
                        f"Stream '{stream_name}' property '{k}' must be numeric."
                    ) from exc
        return stream_ts

    def _run_dynamic_sequential(self, solver, processing_order, results, t_eval, num_steps, start_time, end_time):
        """Standard sequential dynamic simulation for non-recycle flowsheets."""
        # Process units in topological order
        for unit_name in processing_order:
            unit = self.units[unit_name]
            cfg = self.config["units"][unit_name]
            inlet_name = cfg["in"]
            outlets = self._get_unit_outlets(unit_name)
            if isinstance(inlet_name, list):
                inlet_stream_ts = self._merge_inlet_timeseries(inlet_name, results, num_steps, t_eval)
            else:
                inlet_stream_ts = results[inlet_name]

            if hasattr(unit, "run_dynamic"):
                logger.info(f"Simulating dynamic unit {unit_name}")
                sol = unit.run_dynamic(
                    inlet_stream_ts, (start_time, end_time), t_eval, solver
                )
                if len(outlets) == 1:
                    results[outlets[0]] = sol
                else:
                    for out_key, stream_name in zip(["retentate_out", "permeate_out"], outlets):
                        if out_key in sol:
                            results[stream_name] = sol[out_key]
            else:
                logger.debug(f"Processing static unit {unit_name} dynamically")
                sol = self._process_static_unit(
                    unit, inlet_stream_ts, t_eval, num_steps
                )
                if len(outlets) == 1:
                    results[outlets[0]] = sol
                else:
                    for out_key, stream_name in zip(["retentate_out", "permeate_out"], outlets):
                        if out_key in sol:
                            results[stream_name] = sol[out_key]

        self.results = results
        logger.info("Dynamic simulation completed")
        return results

    def _run_dynamic_with_recycle(self, solver, processing_order, results, t_eval, num_steps, start_time, end_time):
        """
        Dynamic simulation with recycle loops using timestep-by-timestep processing.
        
        For each timestep:
        1. Use previous timestep's recycle stream values
        2. Process all units in order for this single timestep
        3. Store results, which become inputs for next timestep
        
        The Tank unit's holdup naturally dampens oscillations and provides stability.
        """
        logger.info("Running dynamic simulation with recycle (timestep-by-timestep)")
        
        # For each unit, we need to track if it has run_dynamic or not
        # For dynamic units like Tank, we need special handling

        # Collect component names up-front so all stream timeseries are initialized consistently.
        components = set()
        for stream_data in self.config["streams"].values():
            z_cfg = stream_data.get("z")
            if isinstance(z_cfg, dict):
                components.update(z_cfg.keys())
        for stream_ts in results.values():
            z_ts = stream_ts.get("z")
            if isinstance(z_ts, dict):
                components.update(z_ts.keys())
        for unit_name in processing_order:
            if self.config["units"][unit_name].get("type") == "Tank":
                unit = self.units[unit_name]
                if hasattr(unit, "initial_n") and isinstance(unit.initial_n, dict):
                    components.update(unit.initial_n.keys())
        
        # Initialize all output streams in results
        for unit_name in processing_order:
            for outlet_name in self._get_unit_outlets(unit_name):
                if outlet_name not in results:
                    # Initialize with zeros/empty structure
                    results[outlet_name] = {
                        "time": t_eval.tolist(),
                        "T": [_DEFAULT_T] * num_steps,
                        "P": [_DEFAULT_P] * num_steps,
                        "flowrate": [0.0] * num_steps,
                        "z": {comp: [0.0] * num_steps for comp in components},
                    }
        
        # Ensure all streams have composition timeseries
        for stream_name in results:
            if "z" not in results[stream_name]:
                results[stream_name]["z"] = {}
            for comp in components:
                if comp not in results[stream_name]["z"]:
                    results[stream_name]["z"][comp] = [0.0] * num_steps
        
        # Store Tank states for dynamic integration
        tank_states = {}
        for unit_name in processing_order:
            if self.config["units"][unit_name]["type"] == "Tank":
                unit = self.units[unit_name]
                tank_states[unit_name] = {
                    "n": {c: unit.initial_n.get(c, 0.0) for c in components},
                    "T": unit.initial_T
                }
        
        # Process timestep by timestep
        for step in range(num_steps):
            t = t_eval[step]
            dt = t_eval[1] - t_eval[0] if num_steps > 1 else 1.0

            if step % 5 == 0:
                logger.debug(f"Processing timestep {step}/{num_steps} (t={t:.2f})")

            # Propagate previous step's tear stream values so upstream units (e.g. tank_1)
            # see the most recently computed recycle values, not the zero initialisation.
            if step > 0:
                for tear_stream in self.tear_streams:
                    if tear_stream in results:
                        prev = step - 1
                        stream_ts = results[tear_stream]
                        for prop, values in stream_ts.items():
                            if prop == "time":
                                continue
                            if prop == "z":
                                for comp_ts in values.values():
                                    comp_ts[step] = comp_ts[prev]
                            elif isinstance(values, list) and step < len(values):
                                values[step] = values[prev]

            # Process each unit for this timestep
            for unit_name in processing_order:
                unit = self.units[unit_name]
                cfg = self.config["units"][unit_name]
                
                # Get merged inlet snapshot for this timestep
                inlet_snapshot = self._get_inlet_snapshot(results, cfg, step)
                
                if cfg["type"] == "Tank":
                    # Special handling for Tank - integrate one step
                    outlet_snapshot, tank_states[unit_name] = self._integrate_tank_step(
                        unit, inlet_snapshot, tank_states[unit_name], dt, components
                    )
                else:
                    # Static unit - just run it
                    outlet_snapshot = unit.run(inlet_snapshot)
                
                # Store results for this timestep
                if "out" in cfg:
                    self._store_snapshot_in_timeseries(results, cfg["out"], outlet_snapshot, step, components)
                else:
                    for out_key in ["retentate_out", "permeate_out"]:
                        if out_key in cfg and out_key in outlet_snapshot:
                            self._store_snapshot_in_timeseries(results, cfg[out_key], outlet_snapshot[out_key], step, components)
        
        self.results = results
        logger.info("Dynamic simulation with recycle completed")
        return results

    @staticmethod
    def _fill_row(
        arr: np.ndarray, idx: int, src_list: list[float], num_steps: int
    ) -> None:
        """Copy up to *num_steps* values from *src_list* into row *idx* of *arr*."""
        n: int = min(len(src_list), num_steps)
        if n:
            arr[idx, :n] = src_list[:n]

    def _merge_inlet_timeseries(
        self,
        inlet_names: list[str],
        results: dict[str, StreamTimeseries],
        num_steps: int,
        t_eval: "np.ndarray",
    ) -> MergedInletTimeseries:
        """Merge multiple inlet timeseries into one by flow-weighted averaging T/P/z and summing flowrate."""
        if not inlet_names:
            raise ValueError("inlet_names must not be empty")

        # Discover all component names across inlets
        all_comps: list[str] = sorted(
            {comp for name in inlet_names for comp in results.get(name, {}).get("z", {})}
        )
        n_inlets = len(inlet_names)

        # Build numpy arrays: shape (n_inlets, num_steps)
        flows = np.zeros((n_inlets, num_steps))
        T_arr = np.full((n_inlets, num_steps), _DEFAULT_T)
        P_arr = np.full((n_inlets, num_steps), _DEFAULT_P)
        z_arrs: dict[str, np.ndarray] = {
            comp: np.zeros((n_inlets, num_steps)) for comp in all_comps
        }

        for idx, name in enumerate(inlet_names):
            ts = results.get(name, {})
            self._fill_row(arr=flows, idx=idx, src_list=ts.get("flowrate", []), num_steps=num_steps)
            self._fill_row(arr=T_arr, idx=idx, src_list=ts.get("T", []), num_steps=num_steps)
            self._fill_row(arr=P_arr, idx=idx, src_list=ts.get("P", []), num_steps=num_steps)
            z_ts = ts.get("z", {})
            for comp in all_comps:
                self._fill_row(arr=z_arrs[comp], idx=idx, src_list=z_ts.get(comp, []), num_steps=num_steps)

        # Vectorized flow-weighted average
        total_flow = flows.sum(axis=0)                       # (num_steps,)
        safe_flow = np.where(total_flow > 0, total_flow, 1.0)

        merged_T = np.where(total_flow > 0, (flows * T_arr).sum(axis=0) / safe_flow, _DEFAULT_T)
        merged_P = np.where(total_flow > 0, (flows * P_arr).sum(axis=0) / safe_flow, _DEFAULT_P)

        merged_z = {}
        for comp in all_comps:
            merged_z[comp] = np.where(
                total_flow > 0, (flows * z_arrs[comp]).sum(axis=0) / safe_flow, 0.0
            ).tolist()

        return MergedInletTimeseries(
            time=t_eval.tolist(),
            T=merged_T.tolist(),
            P=merged_P.tolist(),
            flowrate=total_flow.tolist(),
            z=merged_z,
        )

    def _get_inlet_snapshot(self, results: dict, unit_cfg: dict, step: int) -> dict:
        """Get inlet stream snapshot for a specific timestep, handling multiple inlets."""
        inlets = unit_cfg.get("in")
        if isinstance(inlets, str):
            inlets = [inlets]

        snapshots = [self._extract_snapshot(results.get(name, {}), step) for name in inlets]
        if len(snapshots) == 1:
            return snapshots[0]
        return self._flow_weighted_merge_streams(snapshots)

    def _extract_snapshot(self, stream_ts, step):
        """Extract a single timestep snapshot from a timeseries."""
        snapshot = {}
        for prop, values in stream_ts.items():
            if prop == "time":
                continue
            if prop == "z":
                snapshot["z"] = {
                    comp: comp_ts[step] if step < len(comp_ts) else comp_ts[-1]
                    for comp, comp_ts in values.items()
                }
            else:
                if isinstance(values, list) and step < len(values):
                    snapshot[prop] = values[step]
                elif isinstance(values, list) and values:
                    snapshot[prop] = values[-1]
                else:
                    snapshot[prop] = values
        return snapshot

    def _integrate_tank_step(self, unit, inlet_snapshot, state, dt, components):
        """
        Integrate Tank for one timestep using simple Euler method.
        
        Args:
            unit: Tank unit instance
            inlet_snapshot: Dict with T, P, flowrate, z for inlet
            state: Current Tank state {n: {comp: moles}, T: temperature}
            dt: Timestep size
            components: Set of component names
            
        Returns:
            tuple: (outlet_snapshot, new_state)
        """
        # Current state
        n = state["n"]
        T = state["T"]
        n_tot = max(1e-12, sum(n.values()))
        
        # Inlet properties
        F_in = inlet_snapshot.get("flowrate", 0.0)
        z_in = inlet_snapshot.get("z", {})
        T_in = inlet_snapshot.get("T", _DEFAULT_T)
        
        # Outlet flow (constant from Tank params)
        F_out = unit.outflow
        
        # Current compositions
        x = {c: n.get(c, 0.0) / n_tot for c in components}
        
        # Molar balances: dn/dt = F_in * z_in - F_out * x
        new_n = {}
        for comp in components:
            dn_dt = F_in * z_in.get(comp, 0.0) - F_out * x.get(comp, 0.0)
            new_n[comp] = max(0.0, n.get(comp, 0.0) + dn_dt * dt)
        
        # Simple energy balance (approximate)
        # dT/dt ~ (F_in * T_in - F_out * T) / n_tot (simplified)
        if n_tot > 1e-10:
            dT_dt = (F_in * T_in - F_out * T + unit.duty / 100) / n_tot  # duty scaled
            new_T = T + dT_dt * dt
        else:
            new_T = T_in
        
        new_T = max(200, min(500, new_T))  # Bound temperature
        
        # New state
        new_state = {"n": new_n, "T": new_T}
        
        # Outlet snapshot
        new_n_tot = max(1e-12, sum(new_n.values()))
        outlet_snapshot = {
            "T": new_T,
            "P": unit.P,
            "flowrate": F_out,
            "z": {c: new_n.get(c, 0.0) / new_n_tot for c in components}
        }
        
        return outlet_snapshot, new_state

    def _store_snapshot_in_timeseries(self, results, stream_name, snapshot, step, components):
        """Store a snapshot into the appropriate timestep of a timeseries."""
        if stream_name not in results:
            return
        
        for prop, value in snapshot.items():
            if prop == "z":
                for comp, comp_val in value.items():
                    if comp in results[stream_name]["z"]:
                        results[stream_name]["z"][comp][step] = comp_val
            else:
                if prop in results[stream_name]:
                    results[stream_name][prop][step] = value
