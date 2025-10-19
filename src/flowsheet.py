from loguru import logger

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

        for name, ucfg in self.config["units"].items():
            utype = ucfg["type"]
            cls = unit_types.get(utype)
            if not cls:
                logger.error(f"Unknown unit type: {utype}")
                raise ValueError(f"Unknown unit type: {utype}")
            self.units[name] = cls(name, **{k: v for k, v in ucfg.items() if k != "type"})
            logger.info(f"Built unit {name} of type {utype}")

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
    def run_dynamic(self, solver):
        logger.info("Running dynamic simulation")
        sim = self.config["simulation"]
        duration = int(sim.get("duration", 100))
        timestep = int(sim.get("timestep", 1))
        logger.info(f"Duration: {duration}, Timestep: {timestep}")
        t_eval = list(range(0, duration + 1, timestep))

        # For now assume 1 tank system: feed -> tank -> product
        results = {}
        for stream_name, stream_data in self.config["streams"].items():
            results[stream_name] = stream_data

        for uname, unit in self.units.items():
            cfg = self.config["units"][uname]
            inlet_name = cfg.get("in")
            outlet_name = cfg.get("out")
            inlet = results[inlet_name]

            if cfg["type"] == "Tank":
                sol = unit.run_dynamic(t_eval, [inlet], outlet_flow=inlet["flowrate"])
                results[outlet_name] = {
                    "time": sol.t.tolist(),
                    "level": sol.y[0].tolist(),
                    "unit": "Tank",
                }
                logger.info(f"Simulated dynamic unit {uname}")
            else:
                results[outlet_name] = unit.run(inlet)
                logger.debug(f"Processed static unit {uname}")

        self.results = results
        logger.info("Dynamic simulation completed")
        return results
