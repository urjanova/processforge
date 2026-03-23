# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.14] - 2026-03-23

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
