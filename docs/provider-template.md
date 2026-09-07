# Provider Authoring Guide

This guide shows how to add a new calculation-engine provider to Processforge.
Providers bridge flowsheet stream/unit declarations to a backend library or
service (CoolProp, Cantera, OpenMC, FESTIM, …).

## 1. Choose a base class

### `AbstractProvider` — the minimal contract

Use this when your provider supplies thermodynamic properties and/or unit-level
calculations that plug into the sequential-modular (SM) flowsheet loop.

```python
# processforge/providers/myengine_provider.py
from processforge.providers.base import AbstractProvider
from processforge.providers.registry import register_provider

class MyEngineProvider(AbstractProvider):
    def initialize(self, provider_config, flowsheet_config):
        # Load mechanisms, validate config, open files, etc.
        pass

    def get_thermo_properties(self, stream: dict) -> dict:
        # Return {"H": J/mol, "Cp": J/(mol·K), "K_values": {comp: K}}
        pass

    def compute_unit(self, unit_type: str, config: dict, inlet: dict):
        # Return an outlet stream dict, or None to fall through to default logic.
        return None

    def teardown(self):
        pass

register_provider("myengine", MyEngineProvider)
```

### `BaseSimulationProvider` — for SolverUnit engines

Use this when your provider runs standalone simulations triggered by a
`SolverUnit` (like OpenMC or FESTIM).  It provides default stubs for
`get_thermo_properties`, `compute_unit`, and `teardown`, plus a writable
run-directory helper.

```python
from processforge.providers.base_simulation_provider import BaseSimulationProvider
from processforge.providers.errors import ProviderNotAvailableError, classify_run_error, make_failed_output
from processforge.providers.registry import register_provider
from processforge.types import EngineOutput, OutputProvenance

class MyEngineProvider(BaseSimulationProvider):
    def initialize(self, provider_config, flowsheet_config):
        try:
            import myengine
        except ImportError as exc:
            raise ProviderNotAvailableError("myengine is not installed") from exc
        self._provider_output_dir = self._expand_output_dir(provider_config.output_dir)
        self._materials = dict(flowsheet_config.materials.items())
        self._initialized = True

    def run_simulation(self, unit_config, inlet) -> EngineOutput:
        if not self._initialized:
            raise RuntimeError("MyEngineProvider has not been initialized")
        run_dir = self._resolve_run_dir()
        try:
            result = ...  # run the engine
        except Exception as exc:
            err = classify_run_error("myengine", exc)
            return make_failed_output("myengine", unit_config.sim_type, run_dir, err)
        return EngineOutput(
            status="completed",
            engine="myengine",
            sim_type=unit_config.sim_type,
            fields=[...],
            artifacts=[...],
            diagnostics={"run_dir": str(run_dir)},
            provenance=OutputProvenance(),
            run_dir=str(run_dir),
        )

register_provider("myengine", MyEngineProvider)
```

### `BaseJacobianMixin` — for analytic Jacobian contributions

If your engine can provide cheap analytic or semi-analytic Jacobian blocks for
certain unit types, mix this in alongside `AbstractProvider` or
`BaseSimulationProvider`:

```python
from processforge.providers.base_jacobian_mixin import BaseJacobianMixin

class MyEngineProvider(AbstractProvider, MyEngineJacobianMixin):
    ...
```

`BaseJacobianMixin.compute_jacobian_block` returns `[]` by default, which tells
the caller to fall back to global finite differences for that unit.  Override it
for specific `unit_type` values.  See `CanteraJacobianMixin` for an example.

## 2. Register the provider in the catalog

Add a `ProviderCatalogEntry` to `processforge.providers.registry._PROVIDER_CATALOG`:

```python
"myengine": ProviderCatalogEntry(
    module="processforge.providers.myengine_provider",
    class_name="MyEngineProvider",
    optional_dep="myengine",  # pip extras name, or None
    description="My custom engine provider",
    docker_image=None,        # or a container image for containerized providers
    default_port=None,        # or the container's default API port
),
```

If `optional_dep` is set, users install it with `pip install "processforge[myengine]"`.

## 3. Implement `validate_material` (optional)

If your provider has provider-specific material rules, override the classmethod:

```python
@classmethod
def validate_material(cls, mat_name: str, mat_def, unit_cfg) -> list[str]:
    errors = []
    if mat_def.extra is None or "my_required_field" not in mat_def.extra:
        errors.append(f"Material '{mat_name}' is missing 'extra.my_required_field'.")
    return errors
```

Return an empty list when the material is valid.  These errors are surfaced
by `pf validate` / `pf plan`.

## 4. Add simulation types (SolverUnit engines only)

For engine providers with multiple `sim_type` values, define strategies:

```python
from processforge.providers._sim_strategy import SimStrategy
from processforge.providers.myengine.strategies import (
    MyEngineSimStrategy,
    register_myengine_sim_type,
)

class SteadyState(MyEngineSimStrategy):
    sim_type = "steady_state"

    def build(self, myengine, solver_cfg, geometry_cfg, materials_map, helpers):
        ...
        return model_objects

register_myengine_sim_type("steady_state", SteadyState)
```

Reuse the shared registry via `processforge.providers._sim_strategy` so the
provider dispatch is consistent with OpenMC/FESTIM.

## 5. Containerized vs in-process checklist

| Concern | In-process provider | Containerized provider |
|---|---|---|
| Backend import | In `initialize()` | Only inside the container |
| CLI side | Provider class itself | `ContainerProviderClient` |
| HTTP contract | N/A | `POST /run`, `GET /health` |
| `run_simulation` | Runs engine directly | Runs in container; CLI client receives `EngineOutput` |
| `output_dir` | Resolved locally | Resolved by container server |

For a containerized provider:

1. Implement the provider class (subclassing `BaseSimulationProvider`) inside the
the container image under `processforge.api.serve`.
2. Add `docker_image` and `default_port` to the catalog entry.
3. The CLI side automatically uses `ContainerProviderClient`; no extra code is
needed in `processforge.providers`.

## 6. Error handling

Use the provider exception hierarchy consistently:

* `ProviderNotAvailableError` — backend is missing or container is unreachable.
* `ProviderRuntimeError` — engine crashed during a run.
* `ProviderValidationError` — material/unit config is invalid for this provider.
* `ProviderConfigError` — the `provider_config` block itself is invalid.

For SolverUnit runs, catch engine exceptions and return a classified result:

```python
from processforge.providers.errors import classify_run_error, make_failed_output

except Exception as exc:
    err = classify_run_error("myengine", exc)
    return make_failed_output("myengine", sim_type, run_dir, err)
```

Add engine-specific error signatures to `processforge.providers.errors` if
`myengine` has recognizable failure modes.

## 7. Testing tips

* Subclass tests should check that the provider subclasses `AbstractProvider`
  and implements `initialize`, `teardown`, etc. (see `tests/test_providers.py`).
* Use fake modules / monkeypatching to test strategy `build()` methods without
  installing heavy optional dependencies.
* Test `validate_material` with valid, missing, and malformed inputs.
* If you add error categories, add corresponding tests in
  `tests/test_provider_errors.py`.

## 8. Minimal complete example

```python
"""MyEngine provider — a minimal example."""
from __future__ import annotations

from processforge.providers.base import AbstractProvider
from processforge.providers.errors import ProviderNotAvailableError
from processforge.providers.registry import register_provider


class MyEngineProvider(AbstractProvider):
    def initialize(self, provider_config, flowsheet_config):
        try:
            import myengine
            self._engine = myengine
        except ImportError as exc:
            raise ProviderNotAvailableError(
                "myengine is not installed. Install with: pip install 'processforge[myengine]'"
            ) from exc

    def get_thermo_properties(self, stream: dict) -> dict:
        T, P, z = stream["T"], stream["P"], stream["z"]
        return {
            "H": self._engine.enthalpy(z, T, P),
            "Cp": self._engine.cp(z, T, P),
            "K_values": {c: 1.0 for c in z},
        }

    def compute_unit(self, unit_type: str, config: dict, inlet: dict):
        # MyEngine does not intercept unit calculations.
        return None

    def teardown(self):
        self._engine = None


register_provider("myengine", MyEngineProvider)
```

Add the corresponding config model to `processforge.types` and the provider
catalog entry, and the new backend is available in flowsheets.
