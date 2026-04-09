"""ModelicaProvider — dynamic co-simulation via OpenModelica FMU.

The FMU is a *managed resource*: it is compiled once during ``initialize()``,
loaded as an ``FMU2Slave`` co-simulation object, and stepped at each
``compute_unit()`` call.  ``teardown()`` (or the context manager ``__exit__``)
terminates and frees the FMU instance.

Flowsheet JSON declaration::

    "providers": {
        "openmodelica": {"type": "modelica"}
    }
    "default_provider": "openmodelica"

Install extras: ``pip install "processforge[modelica]"``
                (installs OMPython ≥ 4.0 and fmpy ≥ 0.3)
"""
from __future__ import annotations

import os
import tempfile
from typing import TYPE_CHECKING, Optional

from loguru import logger

from .base import AbstractProvider
from .modelica_jacobian import ModelicaJacobianMixin
from .registry import register_provider

if TYPE_CHECKING:
    from processforge.types import FlowsheetConfig, ModelicaProviderConfig


class ModelicaProvider(AbstractProvider, ModelicaJacobianMixin):
    """Co-simulation provider backed by an OpenModelica-compiled FMU.

    Lifecycle
    ---------
    1. ``initialize()``: transpiles the flowsheet JSON to Modelica (.mo),
       compiles it via OMPython to a Model Exchange / Co-Simulation FMU 2.0,
       then loads it with *fmpy* and enters initialization mode.
    2. ``compute_unit()``: advances the FMU by one ``dt`` step and reads
       output variable values for the requested unit's outlet.
    3. ``teardown()``: terminates the FMU and frees the instance.

    Thermo fall-through
    -------------------
    Modelica does not expose a generic per-stream thermo property service.
    ``get_thermo_properties()`` falls back to CoolProp so that units without
    an explicit ``"provider"`` key still receive correct thermo data.
    """

    def initialize(
        self,
        provider_config: "ModelicaProviderConfig",
        flowsheet_config: "FlowsheetConfig",
    ) -> None:
        config_path = flowsheet_config.extra.get("_config_path")
        if not config_path:
            raise RuntimeError(
                "ModelicaProvider requires '_config_path' in flowsheet_config. "
                "This is set automatically by simulate.py — do not call "
                "ModelicaProvider.initialize() directly without it."
            )

        # Check for required optional packages early.
        try:
            import fmpy  # noqa: F401
            from fmpy import read_model_description, extract
            from fmpy.fmi2 import FMU2Slave
        except ImportError as exc:
            raise RuntimeError(
                "fmpy is not installed. "
                "Install with: pip install 'processforge[modelica]'"
            ) from exc

        from processforge.modelica.transpiler import transpile, _derive_model_name
        from processforge.modelica.omc_runner import compile_modelica

        output_dir = provider_config.output_dir
        os.makedirs(output_dir, exist_ok=True)

        logger.info("ModelicaProvider: transpiling flowsheet to Modelica…")
        mo_path = transpile(config_path, output_dir=output_dir)

        model_name = _derive_model_name(flowsheet_config, config_path)
        logger.info(f"ModelicaProvider: compiling '{model_name}' via OMPython…")
        fmu_path = compile_modelica(mo_path, model_name, output_dir=output_dir)

        logger.info(f"ModelicaProvider: loading FMU '{fmu_path}'…")
        self._model_desc = read_model_description(fmu_path)
        self._fmu_dir = extract(fmu_path)
        self._fmu = FMU2Slave(
            guid=self._model_desc.guid,
            unzipDirectory=self._fmu_dir,
            modelIdentifier=self._model_desc.coSimulation.modelIdentifier,
            instanceName="processforge_cosim",
        )
        self._fmu.instantiate()
        self._fmu.setupExperiment(startTime=0.0)
        self._fmu.enterInitializationMode()
        self._fmu.exitInitializationMode()
        self._t = 0.0

        # Build a quick lookup: variable name → value reference
        self._var_refs: dict[str, int] = {
            v.name: v.valueReference
            for v in self._model_desc.modelVariables
        }
        logger.info("ModelicaProvider: FMU ready for co-simulation.")
        # Attempt to set up a second ME FMU instance for directional derivatives.
        self._initialize_jacobian_fmu(fmu_path, self._model_desc)

    def get_thermo_properties(self, stream: dict) -> dict:
        """Fall back to CoolProp — Modelica does not expose thermo services."""
        from processforge.thermo import get_enthalpy_molar, get_Cp_molar, get_K_values

        z, T, P = stream["z"], stream["T"], stream["P"]
        return {
            "H": get_enthalpy_molar(z, T, P),
            "Cp": get_Cp_molar(z, T, P),
            "K_values": get_K_values(list(z.keys()), T, P),
        }

    def compute_unit(
        self,
        unit_type: str,
        config: dict,
        inlet: dict,
    ) -> Optional[dict]:
        """Advance the FMU one step and extract the unit's outlet stream."""
        if not hasattr(self, "_fmu") or self._fmu is None:
            return None

        dt = config.get("dt", 1.0)
        try:
            self._fmu.doStep(
                currentCommunicationPoint=self._t,
                communicationStepSize=dt,
            )
            self._t += dt
        except Exception as exc:
            logger.warning(f"ModelicaProvider: FMU doStep failed ({exc}). Returning inlet.")
            return dict(inlet)

        # Read outlet variables from FMU — uses Modelica naming convention
        # stream_{T|P|flowrate} for the unit's output stream.
        outlet_stream = config.get("out", "")
        outlet = dict(inlet)
        for prop, key in [("T", f"{outlet_stream}_T"), ("P", f"{outlet_stream}_P"),
                          ("flowrate", f"{outlet_stream}_flowrate")]:
            if key in self._var_refs:
                try:
                    val = self._fmu.getReal([self._var_refs[key]])[0]
                    outlet[prop] = val
                except Exception:
                    pass  # keep inlet value if read fails

        return outlet

    def teardown(self) -> None:
        self._teardown_jacobian_fmu()
        if hasattr(self, "_fmu") and self._fmu is not None:
            try:
                self._fmu.terminate()
                self._fmu.freeInstance()
            except Exception as exc:
                logger.warning(f"ModelicaProvider: FMU teardown error ({exc})")
            finally:
                self._fmu = None
        logger.info("ModelicaProvider: FMU resources released.")


register_provider("modelica", ModelicaProvider)
