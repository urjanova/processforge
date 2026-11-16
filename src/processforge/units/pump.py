import copy

from ..eo.units.pump_eo import PumpEOMixin
from .base import BaseUnitMixin
from .hydraulics import (
    normalize_efficiency,
    pump_power_and_temperature_rise,
    require_pump_properties,
    resolve_properties,
)
from .registry import register_unit


class Pump(BaseUnitMixin, PumpEOMixin):
    """
    Simple steady-state pump model.
    Adds a fixed pressure rise ΔP to the inlet stream and computes the shaft
    power and the adiabatic-inefficiency temperature rise.

    All fluid properties come from the attached thermodynamic provider, or
    from explicit unit parameters when the provider cannot supply them:
    ``cp`` [J/mol/K] and ``molar_volume`` [m³/mol] (or ``density`` [kg/m³]
    together with ``molar_mass`` [kg/mol]). There are no hard-coded fluid
    defaults.

    ``efficiency`` is a fraction in (0, 1]; values greater than 1 are
    interpreted as percent (e.g. ``36`` -> ``0.36``) with a warning.
    """

    def __init__(self, name, delta_p=1e5, efficiency=0.8, **kwargs):
        self.name = name
        self.delta_p = delta_p
        self.efficiency = efficiency
        self.params = {"delta_p": delta_p, "efficiency": efficiency, **kwargs}

    def _run_impl(self, inlet):
        """
        Processes the inlet stream through the pump.

        Outlet pressure is the inlet pressure plus ``delta_p``. Shaft power
        is ``F_mol · v_molar · ΔP / η`` and the lost work heats the fluid:
        ``dT = power · (1 − η) / (F_mol · Cp_molar)``. Molar volume and Cp are
        resolved from the provider or explicit unit parameters (see
        :func:`hydraulics.resolve_properties`).

        Args:
            inlet (dict): Inlet stream with T, P, z, flowrate

        Returns:
            dict: Outlet stream with updated 'P', 'T', 'power', and 'unit'.
        """
        outlet = copy.deepcopy(inlet)
        outlet["P"] = inlet.get("P", 101325.0) + self.delta_p

        eff = normalize_efficiency(self.efficiency, self.name)
        flow_mol = float(inlet.get("flowrate", 1.0))

        provider = getattr(self, "_provider", None)
        props = resolve_properties(
            self.name,
            {
                "z": inlet.get("z", {}),
                "T": inlet.get("T"),
                "P": inlet.get("P", 101325.0),
            },
            self.params,
            getattr(provider, "get_thermo_properties", None) if provider is not None else None,
        )
        require_pump_properties(self.name, props)

        power, dT = pump_power_and_temperature_rise(
            flow_mol, props["molar_volume"], props["cp"], self.delta_p, eff
        )
        outlet["T"] = inlet["T"] + dT

        outlet["power"] = power
        outlet["unit"] = "Pump"
        return outlet


register_unit("Pump", Pump)
