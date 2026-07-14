# Processforge

![processforge-logo](images/processforge-logo.svg)

A lightweight Python framework for process simulation, coupling hydraulic, thermal, and reactor workflows.

## Quick start

```bash
pf init
pf plan flowsheets/hydraulic-chain.json
pf apply flowsheets/hydraulic-chain.json
```

## Install

```bash
uv add processforge
```

For CoolProp-backed units:

```bash
uv add "processforge[coolprop]"
```

Optional solver backends:

```bash
uv add "processforge[eo]"
uv add "processforge[eo-casadi]"
uv add "processforge[modelica]"
```

## Run

Use `plan` to validate and preview changes, then `apply` to solve and store a snapshot.

```bash
pf plan flowsheets/my-flowsheet.json
pf apply flowsheets/my-flowsheet.json
```

For direct run mode:

```bash
pf run flowsheets/my-flowsheet.json
```

## Want more detail?

- `docs/usage.md` — CLI commands, workflows, Docker, and cloud notes
- `docs/flowsheets.md` — flowsheet JSON format, materials, units, and recycle rules
- `flowsheets/` — example flowsheets shipped with the repo

## Python API

```python
from processforge import EOFlowsheet, validate_flowsheet
config = validate_flowsheet("flowsheets/my-flowsheet.json")
fs = EOFlowsheet(config, backend="scipy")
results = fs.run()
```

## License

BSD 3-Clause License
