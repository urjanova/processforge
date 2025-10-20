from loguru import logger
import numpy as np  # Add this import at the top of the file

from .units.pump import Pump
from .units.valve import Valve
from .units.strainer import Strainer
from .units.tank import Tank  # optional
from .solver import Solver


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

    def build_units(self):
        logger.info("Building units")
        unit_types = {
            "Pump": Pump,
            "Valve": Valve,
            "Strainer": Strainer,
            "Tank": Tank,  # optional dynamic element
        }
        for unit_name, unit_config in self.config["units"].items():
            unit_type = unit_config["type"]
            unit_class = unit_types.get(unit_type)
            if not unit_class:
                logger.error(f"Unknown unit type: {unit_type}")
                raise ValueError(f"Unknown unit type: {unit_type}")
            self.units[unit_name] = unit_class(
                unit_name,
                **{
                    k: v
                    for k, v in unit_config.items()
                    if k not in ["type", "in", "out"]
                },
            )
            logger.info(f"Built unit {unit_name} of type {unit_type}")

    def run(self):
        logger.info("Starting simulation")
        self.build_units()
        mode = self.config.get("simulation", {}).get("mode", "steady")
        logger.info(f"Simulation mode: {mode}")
        solver = Solver()
        return self.run_dynamic(solver) if mode == "dynamic" else self.run_steady()

    def run_steady(self):
        logger.info("Running steady-state simulation")
        results = {k: v for k, v in self.config["streams"].items()}
        for uname, unit in self.units.items():
            cfg = self.config["units"][uname]
            inlet = results[cfg["in"]]
            results[cfg["out"]] = unit.run(inlet)
            logger.debug(f"Processed unit {uname}")
        self.results = results
        logger.info("Steady-state simulation completed")
        return results

    def _get_processing_order(self):
        """Determines the processing order of units using a topological sort."""
        # Create a map of stream -> producing unit
        stream_producers = {
            unit_config["out"]: unit_name
            for unit_name, unit_config in self.config["units"].items()
        }

        # Create adjacency list (unit -> downstream units)
        adj = {unit_name: [] for unit_name in self.units}
        in_degree = {unit_name: 0 for unit_name in self.units}

        for unit_name, unit_config in self.config["units"].items():
            inlet_stream = unit_config.get("in")
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
            logger.error(
                "Cycle detected in flowsheet graph. Cannot determine processing order."
            )
            raise ValueError("Cycle detected in the flowsheet graph.")

        logger.info(f"Processing order: {order}")
        return order

    def run_dynamic(self, solver):
        logger.info("Running dynamic simulation")
        sim_config = self.config.get("simulation", {})
        t_start, t_end, dt = (
            sim_config.get("t0", 0),
            sim_config.get("tf", 100),
            sim_config.get("dt", 1),
        )
        t_eval = np.arange(t_start, t_end + dt, dt)
        num_steps = len(t_eval)

        # Initialize results with timeseries structure for all streams
        results = {}
        for stream_name, stream_data in self.config["streams"].items():
            stream_ts = {"time": t_eval.tolist()}
            for k, v in stream_data.items():
                if k == 'z' and isinstance(v, dict):
                    # Create {'comp': [val, val, ...]}
                    stream_ts[k] = {comp: [comp_val] * num_steps for comp, comp_val in v.items()}
                else:
                    # Create {'prop': [val, val, ...]}
                    stream_ts[k] = [v] * num_steps
            results[stream_name] = stream_ts

        processing_order = self._get_processing_order()

        # Process units in topological order
        for unit_name in processing_order:
            unit = self.units[unit_name]
            cfg = self.config["units"][unit_name]
            inlet_name = cfg["in"]
            outlet_name = cfg["out"]
            inlet_stream_ts = results[inlet_name]

            # Check if the unit has a specialized dynamic method
            if hasattr(unit, "run_dynamic"):
                logger.info(f"Simulating dynamic unit {unit_name}")
                # Dynamic units handle their own timeseries output
                sol = unit.run_dynamic(inlet_stream_ts, (t_start, t_end), t_eval, solver)
                # This part may need adjustment based on what `run_dynamic` returns
               # Assuming it returns a full timeseries dictionary for the outlet stream
                results[outlet_name] = sol
            else:
                logger.debug(f"Processing static unit {unit_name} dynamically")        
                # For static units, apply their logic at each time step
                
                # Initialize outlet timeseries structure
                outlet_stream_ts = {prop: [] for prop in inlet_stream_ts if prop != 'z'}
                outlet_stream_ts["time"] = t_eval.tolist()
                
                # Specifically initialize 'z' as a dictionary of lists
                if 'z' in inlet_stream_ts:
                    # Get component names from the keys of the 'z' dictionary
                    outlet_stream_ts['z'] = {comp: [] for comp in inlet_stream_ts['z'].keys()}

                for i in range(num_steps):
                    # Create a snapshot of the inlet stream at time t
                    inlet_snapshot = {}
                    for prop, values in inlet_stream_ts.items():
                        if prop == 'time':
                            continue
                        if prop == 'z':
                            inlet_snapshot['z'] = {comp: comp_ts[i] for comp, comp_ts in values.items()}
                        else:
                            inlet_snapshot[prop] = values[i]
                    
                    outlet_snapshot = unit.run(inlet_snapshot)

                    # Append results to the outlet timeseries
                    for prop, value in outlet_snapshot.items():
                        if prop == 'z':
                            for comp, comp_val in value.items():
                                # Ensure component exists in the structure before appending
                                if comp not in outlet_stream_ts['z']:
                                    outlet_stream_ts['z'][comp] = [0.0] * i # pad with zeros
                                outlet_stream_ts['z'][comp].append(comp_val)
                        else:
                            if prop not in outlet_stream_ts:
                                outlet_stream_ts[prop] = []
                            outlet_stream_ts[prop].append(value)
                
                # Ensure all properties are fully populated. If a static unit didn't
                # modify a property, carry over the value from the inlet.
                for prop, values in outlet_stream_ts.items():
                    if prop == 'time':
                        continue
                    if prop == 'z':
                        for comp, comp_values in values.items():
                            if len(comp_values) < num_steps:
                                # Pad with inlet values if missing
                                inlet_comp_ts = inlet_stream_ts.get('z', {}).get(comp, [0.0] * num_steps)
                                outlet_stream_ts['z'][comp] = inlet_comp_ts
                    elif len(values) < num_steps:
                        # Pad with inlet values if missing
                        outlet_stream_ts[prop] = inlet_stream_ts.get(prop, [0.0] * num_steps)


                results[outlet_name] = outlet_stream_ts

        self.results = results
        logger.info("Dynamic simulation completed")
        return results
