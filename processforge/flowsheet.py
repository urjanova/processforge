from .units.heater import Heater
from .units.flash import Flash

UNIT_MAP = {
    "Heater": Heater,
    "Flash": Flash
}

class Flowsheet:
    def __init__(self, data):
        self.streams = data["streams"]
        self.units = []
        for name, params in data["units"].items():
            utype = params["type"]
            if utype in UNIT_MAP:
                self.units.append(UNIT_MAP[utype](name, params))
            else:
                raise ValueError(f"Unknown unit type: {utype}")

    def run(self):
        streams = dict(self.streams)
        for unit in self.units:
            streams = unit.run(streams)
        return streams
