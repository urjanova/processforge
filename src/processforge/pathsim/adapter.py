"""First-class PathSim adapter for Processforge exported FMUs.

This module is intended to be used with the optional ``pathsim`` dependency::

    pip install "processforge[pathsim]"

Example::

    from processforge.pathsim import ProcessForgeFMU

    plant = ProcessForgeFMU("msre_coupled.fmu")
    controller = pathsim.blocks.PID(...)

    sim = pathsim.Simulation(
        blocks=[controller, plant.block],
        connections=[
            pathsim.Connection(controller.output, plant["openmc_reactor.temperature_default"]),
            pathsim.Connection(plant["openmc_reactor.power"], controller.feedback),
        ],
    )
"""
from __future__ import annotations

import json
import os
from typing import Any


def _default_manifest_path(fmu_path: str) -> str:
    """Return the manifest path next to the FMU, if it exists."""
    base, _ = os.path.splitext(fmu_path)
    candidate = f"{base}_fmu_interface.json"
    return candidate if os.path.exists(candidate) else ""


def _load_manifest(fmu_path: str, manifest_path: str | None) -> dict:
    """Load the ``fmu_interface.json`` manifest.

    If *manifest_path* is ``None``, look for a manifest next to the FMU.
    If no manifest is found, build a minimal manifest by inspecting the FMU
    via ``fmpy`` (if installed).
    """
    if manifest_path is None:
        manifest_path = _default_manifest_path(fmu_path)

    if manifest_path:
        with open(manifest_path, "r", encoding="utf-8") as f:
            return json.load(f)

    # Fallback: parse modelDescription.xml with fmpy.
    try:
        from fmpy import read_model_description
    except ImportError as exc:
        raise RuntimeError(
            "No fmu_interface.json manifest found and fmpy is not installed. "
            "Pass manifest_path explicitly or install processforge[modelica]."
        ) from exc

    md = read_model_description(fmu_path)
    variables = []
    for sv in md.modelVariables:
        variables.append(
            {
                "attr_name": sv.name,
                "logical_name": sv.name,
                "causality": sv.causality or "unknown",
                "variability": sv.variability or "continuous",
                "initial_value": float(sv.start) if sv.start is not None else 0.0,
                "description": sv.description or "",
                "unit": sv.unit or "",
            }
        )
    return {"version": "1.0", "flowsheet_name": md.modelName, "variables": variables}


class ProcessForgeFMU:
    """Logical wrapper around a PathSim FMU block.

    Reads the ``fmu_interface.json`` manifest produced by ``pf export-fmu --pathsim``
    and exposes ports by their logical names (``"openmc_reactor.k_eff"``) instead
    of raw FMI variable names (``"out_openmc_reactor_k_eff"``).
    """

    def __init__(
        self,
        fmu_path: str,
        manifest_path: str | None = None,
        fmu_type: str = "CoSimulationFMU",
        dt: float | None = None,
    ) -> None:
        """Load the FMU and manifest.

        Args:
            fmu_path: Path to the ``.fmu`` file.
            manifest_path: Path to the ``fmu_interface.json`` manifest. If ``None``,
                a manifest next to the FMU is used.
            fmu_type: PathSim block class to use (``"CoSimulationFMU"`` or
                ``"ModelExchangeFMU"``).
            dt: Communication step size for co-simulation. If ``None``, the FMU's
                default experiment step size is used.
        """
        import pathsim.blocks

        self.fmu_path = fmu_path
        self.manifest = _load_manifest(fmu_path, manifest_path)

        start_values: dict[str, float] = {}
        self._by_logical: dict[str, dict] = {}
        self._by_attr: dict[str, dict] = {}
        for var in self.manifest.get("variables", []):
            self._by_logical[var["logical_name"]] = var
            self._by_attr[var["attr_name"]] = var
            if var["causality"] in ("input", "parameter"):
                start_values[var["attr_name"]] = float(var.get("initial_value", 0.0))

        fmu_cls = getattr(pathsim.blocks, fmu_type)
        self.block = fmu_cls(fmu_path, start_values=start_values, dt=dt)

    @property
    def inputs(self) -> list[str]:
        """Logical names of all input ports."""
        return [
            name
            for name, var in self._by_logical.items()
            if var["causality"] == "input"
        ]

    @property
    def outputs(self) -> list[str]:
        """Logical names of all output ports."""
        return [
            name
            for name, var in self._by_logical.items()
            if var["causality"] == "output"
        ]

    @property
    def parameters(self) -> list[str]:
        """Logical names of all parameter ports."""
        return [
            name
            for name, var in self._by_logical.items()
            if var["causality"] == "parameter"
        ]

    def __getitem__(self, logical_name: str) -> Any:
        """Return a PathSim ``PortReference`` for a logical port name.

        The returned object can be passed directly to ``pathsim.Connection``.
        """
        from pathsim.utils.portreference import PortReference

        if logical_name not in self._by_logical:
            raise KeyError(
                f"No FMU port named '{logical_name}'. "
                f"Available: {sorted(self._by_logical)}"
            )
        var = self._by_logical[logical_name]
        attr = var["attr_name"]
        causality = var["causality"]
        if causality not in ("input", "output", "parameter"):
            raise KeyError(f"Cannot connect variable with causality '{causality}'")
        return PortReference(self.block, ports=[attr])

    def value(self, logical_name: str) -> float:
        """Read the current value of an output port by logical name."""
        if logical_name not in self._by_logical:
            raise KeyError(
                f"No FMU port named '{logical_name}'. "
                f"Available: {sorted(self._by_logical)}"
            )
        var = self._by_logical[logical_name]
        if var["causality"] != "output":
            raise ValueError(f"'{logical_name}' is not an output")
        return float(self.block.outputs[var["attr_name"]])

    def __contains__(self, logical_name: str) -> bool:
        return logical_name in self._by_logical

    def __repr__(self) -> str:
        return (
            f"<{type(self).__name__} fmu={self.fmu_path!r} "
            f"inputs={len(self.inputs)} outputs={len(self.outputs)}>"
        )
