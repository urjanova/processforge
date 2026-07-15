# Flowsheet Configuration

Processforge flowsheets are defined as JSON files. Use `simulation.mode` to select steady-state or dynamic behavior.

| `mode` | Solver | Use case |
|--------|--------|----------|
| `"steady"` (default) | EO — global Newton-Raphson | Steady-state without Tank units |
| `"dynamic"` | SM — ODE time-marching | Flowsheets containing Tank units |

```json
{
  "metadata": { "name": "My Flowsheet", "version": "2.0" },
  "materials": {
    "Water":   { "id": 1 },
    "Toluene": { "id": 2 },
    "Steel":   { "id": 3 }
  },
  "material_mixes": {
    "Water_Toluene_Mix": {
      "id": 1,
      "percent_type": "ao",
      "components": [
        { "name": "Water",   "fraction": 0.8 },
        { "name": "Toluene", "fraction": 0.2 }
      ]
    }
  },
  "streams": {
    "feed": { "T": 298.15, "P": 101325, "flowrate": 1.0, "material_mix": 1 }
  },
  "units": {
    "pump_1":  { "type": "Pump",  "in": "feed",       "out": "after_pump", "deltaP": 200000, "efficiency": 0.75, "material": 3 },
    "valve_1": { "type": "Valve", "in": "after_pump",  "out": "product",   "pressure_ratio": 0.5,                "material": 3 }
  },
  "simulation": {
    "mode": "steady",
    "backend": "scipy"
  }
}
```

The optional `backend` key selects the EO solver backend: `"scipy"`, `"pyomo"`, or `"casadi"`. Defaults to `"scipy"`.

## Materials and composition

Stream composition can be defined in two ways:

**Explicit `z` dict** — list component mole fractions directly on the stream:

```json
"feed": { "T": 298.15, "P": 101325, "flowrate": 1.0, "z": { "Water": 0.8, "Toluene": 0.2 } }
```

**`material_mix` reference** — define a reusable mix in the top-level `material_mixes` section and reference it by `id`.

```json
"material_mixes": {
  "Water_Toluene_Mix": {
    "id": 1,
    "percent_type": "ao",
    "components": [
      { "name": "Water",   "fraction": 0.8 },
      { "name": "Toluene", "fraction": 0.2 }
    ]
  }
},
"streams": {
  "feed": { "T": 298.15, "P": 101325, "flowrate": 1.0, "material_mix": 1 }
}
```

Rules:
- `z` and `material_mix` are mutually exclusive on a single stream.
- `id` values must be unique across all mixes.
- Each component `name` in a mix must match a key in the top-level `materials` section.
- When all component fractions are provided they must sum to 1.0.

Every unit also requires a `material` integer field pointing to an `id` in the `materials` section. This identifies the structural material the unit is made of, separate from the fluid composition.

## Available Unit Operations

| Unit Type | Mode | Description | Key Parameters |
|-----------|------|-------------|----------------|
| **Pump** | Steady-state (EO) | Adds pressure rise with efficiency losses | `deltaP`, `efficiency` |
| **Valve** | Steady-state (EO) | Isenthalpic pressure reduction | `pressure_ratio` |
| **Strainer** | Steady-state (EO) | Fixed pressure drop element | `deltaP` |
| **Pipes** | Steady-state (EO) & Dynamic | Laminar flow with friction losses | `delta_p`, `diameter` |
| **Tank** | Dynamic only | Well-mixed molar tank (ODE) | `outlet_flow`, `initial_n`, `initial_T`, `P`, `duty` |
| **Flash** | Steady-state (EO, SciPy backend) | Isothermal flash separator | `P` |
| **Heater** | Steady-state (EO, SciPy backend) | Temperature control unit | `duty`, `flowrate` |

> Flash and Heater use CoolProp internally and are supported on the SciPy backend only.
> Pump, Valve, Strainer, and Pipes support all three backends (SciPy, Pyomo, CasADi).

## Recycle streams

Recycle streams require no special configuration. Any stream produced as the `out` of one unit can be used as the `in` of any other unit — including upstream units.

```json
"units": {
  "tank_1": { "type": "Tank", "in": ["feed", "recycle"], "out": "after_tank", ... },
  "pipe_1": { "in": "after_tank", "out": "recycle", ... }
}
```
