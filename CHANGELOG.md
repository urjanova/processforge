# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.7] - 2026-03-02

### Added

- FMU export capability, the ability to generate a FMU output for a flowsheet.


## [0.2.6] - 2026-03-02

### Added

- A new metadata to the Zarr output to include the version of ProcessForge used for the simulation, aiding in reproducibility and tracking of results.


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
