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
        self.build_units()
        mode = self.config.get("simulation", {}).get("mode", "steady")
        logger.info(f"Simulation mode: {mode}")
        solver = Solver()
        return self.run_dynamic(solver) if mode == "dynamic" else self.run_steady()

    def run_steady(self):
        """
        Runs a steady-state simulation for the flowsheet.
        This method initializes the results dictionary with the configured streams,
        then iterates through each unit in the flowsheet. For each unit, it retrieves
        the inlet stream from the results, runs the unit's simulation, and stores
        the outlet stream back into the results. Logging is performed at the start,
        for each unit processed, and upon completion.
        Returns:
            dict: A dictionary containing the updated stream results after simulation.
        """

        logger.info("Running steady-state simulation")
        results = {k: v for k, v in self.config["streams"].items()}
        for unit_name, unit in self.units.items():
            cfg = self.config["units"][unit_name]
            inlet = results[cfg["in"]]
            results[cfg["out"]] = unit.run(inlet)
            logger.debug(f"Processed unit {unit_name}")
        self.results = results
        logger.info("Steady-state simulation completed")
        return results

    def _get_processing_order(self):
        """
        Determines the processing order of units in the flowsheet using topological sorting.
        This method constructs a directed graph where units are nodes and edges represent
        the flow of streams from producer units to consumer units. It then applies Kahn's
        algorithm (a breadth-first topological sort) to compute a valid processing order.
        If a cycle is detected in the graph (indicating invalid flowsheet configuration),
        it logs an error and raises a ValueError.
        Returns:
            list: A list of unit names in the order they should be processed.
        Raises:
            ValueError: If a cycle is detected in the flowsheet graph, preventing
                determination of a valid processing order.
        """
        
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
        """
        Run a dynamic simulation of the flowsheet over a specified time range.
        This method performs a time-dependent simulation by processing units in topological order.
        For units with a specialized 'run_dynamic' method, it delegates to that method for full
        timeseries handling. For static units, it applies their logic at each discrete time step,
        building timeseries data for inlet and outlet streams.
        Args:
            solver: The solver object (e.g., from scipy.integrate) used for dynamic unit simulations.
                    This is passed to unit-specific 'run_dynamic' methods if available.
        Returns:
            dict: A dictionary containing timeseries results for all streams in the flowsheet.
                  Each stream's data includes a 'time' key with the evaluation times, and other
                  properties as lists or nested dictionaries (e.g., 'z' for compositions).
        Notes:
            - Simulation parameters (t0, tf, dt) are retrieved from self.config["simulation"].
            - Results are stored in self.results and returned.
            - For static units, properties not modified by the unit are carried over from the inlet.
            - Assumes units have a 'run' method for static snapshots and optionally 'run_dynamic' for dynamics.
        """

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
