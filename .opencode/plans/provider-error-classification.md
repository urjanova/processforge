# Plan: Robust provider run-error catching & classification (v0.3.17)

Status: PLAN — edits are blocked while in plan mode. Apply after approval.

## Goal
Catch, classify, and document run-time simulation failures that originate **inside
the provider/engine** (e.g. OpenMC nuclear-data / MPI errors) when running in
Docker, attributing them to `source="provider"` (distinct from flowsheet/setup
errors) and surfacing them loudly to the CLI user.

Decisions (confirmed with user):
- Container keeps returning HTTP 200 + structured `EngineOutput` (classification preserved).
- CLI (`pf run`) detects `status="failed"` and exits non-zero with a clear message.
- Scope: OpenMC + FESTIM (shared utility, reusable by any future provider).

## 1. New file `src/processforge/providers/errors.py`
Shared classification utility. Contents:
- Category string constants: `NUCLEAR_DATA`, `CROSS_SECTIONS`, `MPI_ABORT`,
  `GEOMETRY`, `TALLY`, `CONVERGENCE`, `INPUT_VALIDATION`, `ENVIRONMENT`, `UNKNOWN`.
- `_ERROR_SIGNATURES`: ordered `(category, regex, hint)` list. Key signatures:
  - `Nuclear data library does not contain cross sections` → `NUCLEAR_DATA`
    (hint: lower material temp / enable `Settings.temperature` / multipole).
  - cross-section file/path patterns → `CROSS_SECTIONS`.
  - `MPI_ABORT` → `MPI_ABORT` (root cause is the preceding engine error).
  - geometry: `cannot find cell`, `particle lost`, `universe ... not found`, … → `GEOMETRY`.
  - `tally error` / `invalid filter` → `TALLY`.
  - `did not converge` / `maximum number of iterations` → `CONVERGENCE`.
  - `ValueError`/`KeyError`/`TypeError`/`validation error` → `INPUT_VALIDATION`.
  - `permission denied` / `no such file` / `.so` missing → `ENVIRONMENT`.
- `class ProviderRunError(BaseModel)`: `category`, `source="provider"`, `message`
  (concise one-line), `type` (exception class name), `detail` (full text), `hint`.
  - `from_exception(exc, captured="", category=None)` classmethod.
- `classify_run_error(engine, exc, captured="")` → `ProviderRunError` (prefixes message with `[engine]`).
- `make_failed_output(engine, sim_type, run_dir, err, unit="")` → `EngineOutput`
  with `status="failed"`, `error=err`, and `diagnostics` carrying `run_dir`,
  `error`, `error_category`, `error_source`.
- Helpers `_classify_text`, `_hint_for`, `_summarize` (prefers engine `ERROR:` lines).

`errors.py` must NOT import `processforge.types` at module top (only lazily inside
`make_failed_output`) to avoid an import cycle.

## 2. `src/processforge/types.py` — add `error` field to `EngineOutput`
After line 475 (`provenance: ...`), add:
```python
    from processforge.providers.errors import ProviderRunError
    error: Optional[ProviderRunError] = None
```
(Move the `ProviderRunError` import to the top-level imports if a circular import
surfaces; `errors.py` deliberately avoids importing `types`.)

## 3. `src/processforge/providers/openmc_provider.py`
- Add to the `processforge.types` import block: `ProviderRunError` (or import from
  `processforge.providers.errors`).
- In `run_simulation` catch block (~L805): replace the manual
  `EngineOutput(status="failed", diagnostics={...})` with:
  ```python
  err = classify_run_error("openmc", exc)
  return make_failed_output("openmc", sim_type, run_dir, err)
  ```
- (Optional, recommended) capture `model.run` output: keep using `str(exc)` — openmc's
  `executor._run` already raises `RuntimeError` containing the full stderr banner +
  `MPI_ABORT`, so classification works without extra plumbing.

## 4. `src/processforge/providers/festim_provider.py`
- Same change in the `run_simulation` catch block (~L770):
  ```python
  err = classify_run_error("festim", exc)
  return make_failed_output("festim", sim_type, run_dir, err)
  ```

## 5. `src/processforge/cli/run.py` — loud failure
After `results = fs.run()` (steady branch ~L107 and dynamic branch ~L88), add a
failed-output check:
```python
from processforge.types import EngineOutput
_failed = [
    name for name, out in getattr(fs, "engine_outputs", {}).items()
    if isinstance(out, EngineOutput) and out.status == "failed"
]
if _failed:
    for name in _failed:
        e = fs.engine_outputs[name].error
        logger.error(
            f"Unit '{name}' simulation FAILED "
            f"[{getattr(e, 'category', 'unknown')}]: {getattr(e, 'message', '')}"
        )
        if getattr(e, "hint", ""):
            logger.error(f"  hint: {e.hint}")
    raise SystemExit(1)
```
(Place before `collect_outputs` so a failed run is not recorded as success.)

## 6. Docs — new `docs/provider-errors.md`
Document the error taxonomy: categories, `source="provider"` attribution vs
flowsheet validation, example messages (the Ni58 case), and remediation hints.
Link from `docs/provider-docker-images.md`.

## 7. Tests
- `tests/test_provider_errors.py`: `classify_run_error` returns correct category
  for: Ni58 nuclear-data string → `NUCLEAR_DATA`; `MPI_ABORT` text → `MPI_ABORT`;
  geometry/tally strings; unknown fallback → `UNKNOWN`; `from_exception` sets `source`,
  `type`, `hint`; `make_failed_output` yields `status="failed"` with `error` populated
  and `diagnostics.error_category`.
- Extend `tests/test_openmc_provider.py::test_run_exception_returns_failed_result`:
  set `_CTRL.raise_msg` to the real Ni58 message and assert
  `result.error.category == "nuclear_data"` and `result.error.source == "provider"`.
- Extend FESTIM test similarly if a run-failure test exists.

## 8. Version bump + CHANGELOG
- `pyproject.toml`: `version = "0.3.16"` → `version = "0.3.17"`.
- `CHANGELOG.md`: insert new section above `[0.3.16]`:

```markdown
## [0.3.17] - 2026-08-25

### Added
- **Provider run-error classification** (`processforge.providers.errors`): run-time
  failures raised inside a provider/engine (OpenMC, FESTIM, …) are now caught and
  classified into a typed `ProviderRunError` (`source="provider"`) with a category
  (e.g. `nuclear_data`, `cross_sections`, `mpi_abort`, `geometry`, `tally`,
  `convergence`) plus a concrete remediation hint. The `EngineOutput` returned on
  failure now carries this via an `error` field (and `diagnostics.error_category`).

### Changed
- OpenMC and FESTIM providers return the structured `ProviderRunError` instead of a
  bare `diagnostics["error"]` string when a run fails in the container.
- `pf run` now detects any `SolverUnit` whose run returned `status="failed"` and
  exits non-zero with the categorized error and hint, instead of silently saving a
  successful manifest.

### Docs
- Added `docs/provider-errors.md` documenting the provider run-error taxonomy and
  remediation hints.
```

## Verification
- `pytest tests/test_provider_errors.py tests/test_openmc_provider.py tests/test_festim_provider.py`
- `python -c "import processforge.providers.errors"` sanity import.
- `pf run flowsheets/openmc/msre_eigenvalue.json` against the container →
  categorized `nuclear_data`/`cross_sections` error printed, non-zero exit.
