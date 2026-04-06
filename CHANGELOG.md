# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
## [0.2.20] - 2026-04-06

### Added
- Added the `pf init` CLI command to initialise the `.processforge/` project directory and `outputs/` folder.
- Added versioned snapshots to `.pfstate` State Management. Every successful `pf apply` now saves a new numbered snapshot in `snapshots/` with a `latest` file pointer, allowing for safe rollbacks.
- Added structural diff generation (unit additions, removals, modifications) to the `pf plan` and `pf apply` commands to preview topology changes against the last saved state.
- Added divergence guardrails to `pf apply`: if both standard Newton solve and homotopy fallback fail, it automatically reverts to the last good state snapshot and writes a detailed debug report (`*_divergence.json`).
- Added support in `pf run` (dynamic mode) to automatically load the latest converged `.pfstate` values as the initial conditions for time-integration.

### Changed
- Updated `README.md` to reflect new `pf` CLI commands and added details on the new `plan` and `apply` workflow, including snapshot versioning and divergence guardrails.
- Updated `flowsheet.py` to store the solver convergence status and statistics after solving.
- Refactored `StateManager` to fully support Zarr-backed versioning, structural diff detection (`detect_structural_diff`), rollback capability, and state to stream dict conversion.
- Updated `pf apply` to automatically fall back to a cold start when flowsheet topology changes are detected.

## [0.2.19] - 2026-04-06
### Added
- Added `.pfstate` State Management backed by Zarr to persistently store converged simulation values alongside configuration.
- Added drift detection to intelligently compare the "Desired State" (current flowsheet definition) against the "Actual State" (.pfstate).
- Added the `pf apply` CLI command to execute simulations with automatic drift detection and warm-starting.
- Implemented a Homotopy Solver that acts as a step-wise parameter continuation method (10 steps) if the standard non-linear solver fails.

## [0.2.18] - 2026-04-06
### Added
- Added a provider Jacobian contribution API with `JacobianContributor`, `ReferenceState`, and `ReferenceStateRegistry`.
- Added provider mixins for analytic/semi-analytic Jacobian contributions, including `CanteraJacobianMixin` and `ModelicaJacobianMixin`.

### Changed
- Updated global Jacobian evaluation to route provider-contributed Jacobian blocks through the analytic path and fall back to block-sparse finite differences when needed.

## [0.2.17] - 2026-04-06
### Added
- Added the `plan` CLI command for flowsheet validation, DOF analysis, and Mermaid diagram generation.
- Added `validate_flowsheet_dict` utility for validating flowsheet configurations from dictionaries.
- Added `pcl` optional dependency group containing `pint`.

### Modified
- Replaced the `validate` CLI command with the new `plan` command.

## [0.2.16] - 2026-04-05

### Added

- Added the `pf` command to be used interchangably and updated Readme.
- Changed provider behavior: No providers block means all units silently get `coolprop` (zero change for existing flowsheets).
- Added support for `"default_provider": "cantera"` where all units without an explicit "provider" key will use `cantera`.
- Changed `simulation.backend` default to `"scipy"` (was `"pyomo"` in schema default, now resolves to a scipy fallback in code).
- Updated dependencies: `ompython` moved from core dependencies to the `[modelica]` optional group; `fmpy` and `cantera` were added as new optional dependency groups.

## [0.2.15] - 2026-03-23

### Added

- Added validation for `friendly_material_mix_id` in flowsheet JSON, ensuring it is required and unique across all material mixes. This enhances data integrity and consistency in flowsheet definitions.
- Removed the `z` key from the `feed` stream in `closed-loop-chain.json`, as the composition is now defined via `material_mix: 1`. This eliminates redundancy and potential confusion in stream definitions.

## [0.2.12] - 2026-03-16

### Fixed

- Made material mandatory for all unit types. 

## [0.2.11] - 2026-03-16

### Fixed

- Updated `build_units` function to accomodate materials and material_mixes so that it can be processed. 

## [0.2.10] - 2026-03-16

### Added

- Added support for OpenMC material definitions, allowing users to define materials using OpenMC's format for use in simulations. We dont use OpenMC in the flow but this allows users to easily convert their OpenMC material definitions for use in ProcessForge, enhancing compatibility and ease of use for those familiar with OpenMC.

## [0.2.7] - 2026-03-02

### Added

- FMU export capability, the ability to generate a FMU output for a flowsheet.


## [0.2.6] - 2026-03-02

### Added

- A new metadata to the Zarr output to include the version of ProcessForge used for the simulation, aiding in reproducibility and tracking of results.



## [0.2.5] - 2026-03-01

### Removed

- Removed validation to a new [processforge-validation](https://www.github.com/urjanova/processforge-validation) repository, updated dependencies

## [0.2.4] - 2026-02-25

### Fixed

- Fixed a bug in the pump.py file when flowrate is 0.0, both power and mass_flow are zero, making dT = 0/0. Physically, zero flow means no heat exchange, so dT should be 0.0.

## [0.2.3] - 2026-02-25

### Removed

- Removed Upload logs to S3 and upload Zarr to S3 features to simplify codebase and focus on core simulation functionality.

## [0.2.2] - 2026-02-25

### Added

- Upload logs to S3 

## [0.2.1] - 2026-02-25

### Added

- Upload to S3 and set environment variables

## [0.2.0] - 2026-02-25

### Added

- Zarr support for storing simulation results, enabling efficient access and analysis of large datasets.
- Added support for CasADi backend for solving dynamic optimization problems, providing an alternative to Pyomo + IPOPT.

### Fixed
- Fixed a bug in the flowsheet diagram generation where certain unit types were not rendered correctly.

### Changed

- Made image generation optional to speed up simulations and reduce output size when not needed.
- Updated documentation to reflect changes in command-line interface and new features.

### Removed

- CSV / JSON output formats to simplify codebase and focus on Zarr for results storage.


[0.2.0]: https://github.com/urjanova/processforge/compare/v0.1.0...HEAD
