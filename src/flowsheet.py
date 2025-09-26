from .units.heater import Heater
from .units.flash import Flash
from .units.tank import Tank

UNIT_MAP = {
    "Heater": Heater,
    "Flash": Flash,
    "Tank": Tank
}

class Flowsheet:
    def __init__(self, data):
        # top-level streams and units
        self.streams = data.get("streams", {})
        self.unit_defs = data.get("units", {})
        self.simulation = data.get("simulation", {"mode": "steady"})
        # instantiate units
        self.units = []
        for name, params in self.unit_defs.items():
            utype = params["type"]
            if utype in UNIT_MAP:
                self.units.append( (name, UNIT_MAP[utype](name, params)) )
            else:
                raise ValueError(f"Unknown unit type: {utype}")

    def run_steady(self):
        streams = dict(self.streams)
        # run units in declared order (MVP simplification)
        for name, unit in self.units:
            # units operate on streams in-place or return updated dict
            streams = unit.run(streams)
        return streams

    def run_dynamic(self):
        """
        Very simple dynamic orchestration:
        - Expects simulation dict with t0, tf, dt
        - Only integrates Tank units using their inlet streams from self.streams.
        - For each time point we currently only return the tank state snapshots (other
          units could be evaluated per time step in a more complete implementation).
        """
        sim = self.simulation
        t0 = sim.get("t0", 0.0)
        tf = sim.get("tf", 100.0)
        dt = sim.get("dt", 1.0)
        t_eval = [t0 + i*dt for i in range(int((tf - t0)/dt) + 1)]

        # Find tanks
        tanks = [(name, u) for name, u in self.units if isinstance(u, Tank)]
        if not tanks:
            raise ValueError("No Tank units found for dynamic simulation")

        # For MVP: support single tank integration and return its time series.
        # If multiple tanks or interactions are needed, orchestration must be extended.
        name, tank = tanks[0]
        # pass the current flowsheet streams (inlet streams should exist)
        timeseries = tank.run_dynamic(self.streams, (t0, tf), t_eval)

        # Optionally: evaluate steady units each time using current tank state (not done here)
        return timeseries

    def run(self):
        mode = self.simulation.get("mode", "steady")
        if mode == "steady":
            return self.run_steady()
        elif mode == "dynamic":
            return self.run_dynamic()
        else:
            raise ValueError(f"Unknown simulation mode: {mode}")
