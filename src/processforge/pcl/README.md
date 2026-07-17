# PCL — Process Configuration Language

PCL is a TOML-based DSL for defining processforge flowsheets. Files use the `.pcl` extension.

## Syntax

A PCL file has two top-level sections:

```toml
[units]
[streams]
```

### Units

Units reference a process type via the `source` key using the `pf.*.*` namespace:

```toml
[units.my_pump]
source = "pf.hydraulics.pump"
inlet = "feed"
outlet = "discharge"

[units.my_tank]
source = "pf.storage.tank"
inlet = "discharge"
outlet = "product"
```

### Streams

Streams carry material between units:

```toml
[streams.feed]
inlet = { unit = "source", port = "out" }
outlet = { unit = "my_pump", port = "in" }
units_annotation = "kg/s"

[streams.discharge]
inlet = { unit = "my_pump", port = "out" }
outlet = { unit = "my_tank", port = "in" }
units_annotation = "kg/s"
```

## Compilation

Load and compile a `.pcl` file with:

```python
from processforge.pcl import load_pcl, compile_to_dict, PCLCompileError

config = load_pcl("flowsheets/my_design.pcl")
# Returns a dict compatible with processforge flowsheet config
```

The compiler transforms:
- `units[name].source` → `units[name].type` (via the namespace map)
- `streams[name].units_annotation` → `streams[name]._units`

## Namespace

| Source reference | Unit type |
|---|---|
| `pf.hydraulics.pump` | Pump |
| `pf.hydraulics.valve` | Valve |
| `pf.hydraulics.strainer` | Strainer |
| `pf.hydraulics.pipes` | Pipes |
| `pf.storage.tank` | Tank |
| `pf.separations.flash` | Flash |
| `pf.thermal.heater` | Heater |
| `pf.reactors.cstr` | CSTR |
| `pf.reactors.pfr` | PFR |
| `pf.reactors.ideal_gas` | IdealGasReactor |

See `namespace.py` for the authoritative mapping.
