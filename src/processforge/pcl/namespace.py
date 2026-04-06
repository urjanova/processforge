"""Mapping from pf.* namespace references to unit type strings."""

PCL_NAMESPACE: dict[str, str] = {
    "pf.hydraulics.pump":     "Pump",
    "pf.hydraulics.valve":    "Valve",
    "pf.hydraulics.strainer": "Strainer",
    "pf.hydraulics.pipes":    "Pipes",
    "pf.storage.tank":        "Tank",
    "pf.separations.flash":   "Flash",
    "pf.thermal.heater":      "Heater",
    "pf.reactors.cstr":       "CSTR",
    "pf.reactors.pfr":        "PFR",
    "pf.reactors.ideal_gas":  "IdealGasReactor",
}
