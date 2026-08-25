# Provider run errors

When a simulation is executed **inside** a provider/engine (OpenMC, FESTIM, … —
typically inside a Docker container), the run can fail for reasons that live
entirely in the engine. These are different from *flowsheet/setup* errors
(bad JSON, unknown material, schema validation), which are caught earlier during
initialization or flowsheet validation.

Processforge classifies provider run-time failures so they can be distinguished,
surfaced, and remediated:

* Every failed run returns an `EngineOutput` with `status="failed"` and a
  populated `error` field (a `ProviderRunError`).
* The `error.source` is always `"provider"`, so the rest of the framework (and
  you) can tell a runtime engine error apart from a flowsheet configuration error.
* The `error.category` names the failure family; `error.message` is a concise
  one-line summary; `error.detail` is the full captured engine output; and
  `error.hint` is a concrete remediation suggestion.

## Error categories

| Category | Typical cause | Hint |
| --- | --- | --- |
| `nuclear_data` | A nuclide is requested at a temperature absent from the cross-section library (e.g. Ni58 at 400 K when only 300 K is tabulated). | Set the material temperature to an available value, or enable `openmc.Settings.temperature` handling (temperature_method / multipole). |
| `cross_sections` | The `cross_sections` path is wrong/unmounted or the XML is unreadable. | Verify the provider `cross_sections` path and that the data directory is mounted in the container. |
| `mpi_abort` | The solver process crashed (`MPI_ABORT`). | Usually a symptom of a preceding nuclear-data/geometry error — fix that first. Otherwise check container memory/MPI. |
| `geometry` | Invalid geometry, or the source/particles fall outside it. | Check `geometry_config` dimensions and `source_point`, and that every region is filled. |
| `tally` | A tally or filter is invalid / references a missing cell or mesh. | Check `mesh_tallies` / tally scores and filter IDs. |
| `convergence` | The solve did not converge. | Increase batches/iterations or relax tolerances. |
| `input_validation` | The engine rejected a resolved config value (out of range / mistyped). | Inspect `solver_config` / `geometry_config`. |
| `environment` | Permissions, missing shared library, or unwritable output dir. | Check container mounts and output-dir writability. |
| `unknown` | No signature matched. | Inspect the engine log / run directory. |

## Example

The MSRE eigenvalue flowsheet declares `inor` at 400 K, but the FENDL library
ships Ni58 at 300 K. OpenMC aborts inside the container:

```
ERROR: Nuclear data library does not contain cross sections for Ni58 at or near
       400.000000 K. Available temperatures are 300 K.
```

Classified as `nuclear_data`, `pf run` prints:

```
ERROR  Unit 'openmc_solver' simulation FAILED [nuclear_data]: [...] Nuclear data
       library does not contain cross sections for Ni58 ...
ERROR    hint: The cross-section library lacks data at the requested temperature.
       Either set the material temperature to a value present in the library
       (e.g. 300 K), or enable openmc.Settings.temperature handling ...
```

and exits non-zero — the failure is no longer silently recorded as success.

## Implementation

`processforge.providers.errors` provides:

* `ProviderRunError` — the typed record above.
* `classify_run_error(engine, exc, captured="")` — turn an exception (+ optional
  captured stdout/stderr) into a `ProviderRunError`.
* `make_failed_output(engine, sim_type, run_dir, err)` — build the
  `EngineOutput(status="failed")` providers return.

Providers call these in their `run_simulation` `except` block; the CLI
(`pf run`) detects any `engine_outputs` entry with `status="failed"` and exits
non-zero after printing the category and hint.
