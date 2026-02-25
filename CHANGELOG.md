# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
