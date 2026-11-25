# Unit Operations Reference

This page documents every unit type that can appear in a processforge flowsheet's `units` block, with the configurable **properties** — and which properties are computed by the **built-in (SM)** logic versus the **IDAES provider** when a unit declares `"provider": "idaes"`. A new **Hydraulic network mode** is also documented — when enabled in `simulation.hydraulics`, the native SciPy solver determines all stream pressures and element flows simultaneously using physically consistent pressure–flow laws instead of the legacy fixed-value treatment.

## Property vectors

Processforge deliberately keeps **one unit type name with two property vectors**:

- **Built-in (SM)** — the default sequential-modular models in
  `processforge.units`. These units read and compute properties using
  water-like assumptions (ρ = 1000 kg/m³, Cp = 4180 J/kg·K, MW = 0.018 kg/mol)
  unless a thermodynamic provider is attached.
- **IDAES (provider)** — when a unit sets `"provider": "idaes"`, the IDAES
  provider intercepts the same unit types and can compute some of them with
  IDAES steady-state models. IDAES currently supports a core set of properties
  out of the box; more robust property support will be added later.
- **Hydraulic network mode** — opt-in via `simulation.hydraulics.enabled: true`.
  A native nonlinear solver (`scipy.optimize.least_squares`) solves stream
  pressures and element flows each timestep using industry-standard laws:
  IDAES valve law, Hagen–Poiseuille / Darcy–Weisbach pipes, EPANET pump curves,
  and linear design-point resistances for strainers.

Every unit type accepts the same shared keys regardless of vector:

| Key | Type | Description |
|-----|------|-------------|
| `name` | str | Unique unit name. |
| `type` | str | Unit type — see the tables below. |
| `in` | str \| [str] | Inlet stream name(s). |
| `out` | str \| [str] | Outlet stream name(s). |
| `provider` | str | *Optional.* Provider to delegate to (e.g. `"idaes"`). |
| `material` | int | *Optional.* Material `id` the unit is made of. |

Stream state uses: `T` in K, `P` in Pa, `flowrate` in mol/s, `z` a
mole-fraction dict. The built-in models assume water-like properties (ρ = 1000
kg/m³, Cp = 4180 J/kg·K, MW = 0.018 kg/mol) unless a thermodynamic provider is
attached.

## Unit operation properties

The **Built-in (SM)** column indicates whether the default sequential-modular
logic reads the key. The **IDAES (provider)** column indicates whether the
IDAES provider understands the property out of the box. The **Hydraulic**
column indicates whether the key is used by the pressure-driven network solver.

Legend:
- ✅ — property is read on this vector.
- — — property is not modelled on this vector.
- 🔽 — unit type falls through to the built-in logic.
- ⚠️ — legacy property; ignored in hydraulic mode with a warning.

### Pump

Adds a pressure rise to the inlet stream. In hydraulic mode the rise comes
from a `pump_curve`; otherwise a fixed `delta_p` is used.

| Property | Built-in (SM) | IDAES (provider) | Hydraulic | Default | Units | Description |
|----------|---------------|------------------|-----------|---------|-------|-------------|
| `delta_p` | ✅ | ✅ | ✅ (fallback) | `1e5` | Pa | Fixed pressure rise. |
| `efficiency` | ✅ | ✅ | — | `0.8` | — | Adiabatic efficiency. |
| `pump_curve` | — | — | ✅ | — | — | EPANET-style head/ΔP–flow curve (see below). |

**Pump curve format** (`pump_curve`):

```json
{
  "flow": [0.0, 0.0005556, 0.0011112],
  "delta_p": [4000000.0, 3000000.0, 1500000.0]
}
```

or head-based (converted using fluid density from the provider):

```json
{
  "flow": [0.0, 1.0, 2.0],
  "head": [40.8, 30.6, 15.3]
}
```

The curve is fitted to `dP = A − B·F^C` following WNTR conventions:
- 1 point: `A = 4H/3`, `B = H/(3Q²)`, `C = 2`
- 2 points: straight line (`C = 1`)
- 3+ points: `scipy.optimize.curve_fit`

### Valve

Isenthalpic pressure reduction. In **hydraulic mode** the flow is governed by
the IDAES valve law `F² = Cv² · (P_in − P_out) · f(opening)²`. In legacy mode
a fixed `pressure_ratio` is applied.

| Property | Built-in (SM) | IDAES (provider) | Hydraulic | Default | Units | Description |
|----------|---------------|------------------|-----------|---------|-------|-------------|
| `pressure_ratio` | ✅ | ✅ | ⚠️ (ignored) | `0.5` | — | Legacy: P_out = ratio·P_in. |
| `Cv` | — | — | ✅ | — | mol/s/Pa⁰·⁵ | IDAES flow coefficient (hydraulic mode). |
| `valve_opening` | — | — | ✅ | `1.0` | — | Fraction open [0, 1]. |
| `valve_function` | — | — | ✅ | `"linear"` | — | `"linear"`, `"quick_opening"`, `"equal_percentage"`. |
| `alpha` | — | — | ✅ | `100.0` | — | Equal-percentage rangeability. |
| `check_valve` | — | — | ✅ | `false` | — | If `true`, permits only forward flow. |

**Valve functions:**
- `linear`: `f(x) = x`
- `quick_opening`: `f(x) = sqrt(x)`
- `equal_percentage`: `f(x) = alpha^(x−1)`

### Strainer

Fixed pressure drop element. In hydraulic mode it uses a linear design-point
resistance `dP = delta_p_design · (F/F_design)`.

| Property | Built-in (SM) | IDAES (provider) | Hydraulic | Default | Units | Description |
|----------|---------------|------------------|-----------|---------|-------|-------------|
| `delta_p` | ✅ | ✅ | ✅ | `5000.0` | Pa | Pressure drop. |
| `design_flowrate` | — | — | ✅ | — | mol/s | Design flow for resistance scaling. |

### Pipes

Pipe segment with pressure drop. In hydraulic mode it can use full geometry
(Hagen–Poiseuille for laminar, Darcy–Weisbach with Swamee–Jain for turbulent)
or a linear design-point fallback (`delta_p` + `design_flowrate`).

| Property | Built-in (SM) | IDAES (provider) | Hydraulic | Default | Units | Description |
|----------|---------------|------------------|-----------|---------|-------|-------------|
| `delta_p` | ✅ | ✅ | ✅ (fallback) | `1000.0` | Pa | Fixed pressure drop. |
| `diameter` | ✅ | ✅ | ✅ (geometry) | `0.1` | m | Pipe inside diameter. |
| `length` | — | — | ✅ (geometry) | `1.0` | m | Pipe length. |
| `roughness` | — | — | ✅ (geometry) | `0.0` | m | Absolute roughness (`0` = smooth). |
| `design_flowrate` | — | — | ✅ (fallback) | — | mol/s | Design flow for linear scaling. |
| `viscosity` | — | — | ✅ (override) | — | Pa·s | Dynamic viscosity (overrides provider). |
| `molar_volume` | — | — | ✅ (override) | — | m³/mol | Molar volume (overrides provider). |
| `density` | — | — | ✅ (override) | — | kg/m³ | Density (overrides provider). |

### Tank

Well-mixed molar tank. **Dynamic only** (requires dynamic mode).
In hydraulic mode the tank acts as a **pressure boundary**: its inlet and outlet
flows are independent and determined by the network. The inventory is integrated
by the flowsheet driver.

| Property | Built-in (SM) | IDAES (provider) | Hydraulic | Default | Units | Description |
|----------|---------------|------------------|-----------|---------|-------|-------------|
| `outlet_flow` | ✅ | — | ⚠️ (ignored) | `1.0` | mol/s | Legacy outlet flow. |
| `P` | ✅ | — | ✅ (boundary) | `101325.0` | Pa | Operating pressure (fixed). |
| `duty` | ✅ | — | — | `0.0` | W | Heat duty (positive heats). |
| `initial_n` | ✅ | — | — | `{}` | mol | Initial molar inventory. |
| `initial_T` | ✅ | — | — | `300.0` | K | Initial temperature. |
| `flow_in` | ✅ | — | — | `1.0` | mol/s | Inlet flow when stream lacks one. |

### Flash

Isothermal flash separator (Rachford–Rice) producing vapor and liquid
outlets. 🔽 **Fall-through** — needs K-values from a thermo provider.

| Property | Built-in (SM) | IDAES (provider) | Hydraulic | Default | Units | Description |
|----------|---------------|------------------|-----------|---------|-------|-------------|
| `P` | ✅ | — | — | — | Pa | Operating pressure. |
| `in` | ✅ | — | — | — | — | Inlet stream name. |
| `out_vap` | ✅ | — | — | — | — | Outlet vapor stream name. |
| `out_liq` | ✅ | — | — | — | — | Outlet liquid stream name. |

### Heater, CSTR and PFR

These units are **not** supported in IDAES currently. CSTR and PFR require a
Cantera provider (`"provider": "cantera"`). Heater relies on CoolProp for
enthalpy calculations and is not yet IDAES‑compatible.

## Hydraulic network mode

### Enabling

Opt-in via the `simulation` block:

```json
{
  "simulation": {
    "mode": "dynamic",
    "hydraulics": {
      "enabled": true,
      "tol": 1e-8,
      "max_iter": 200,
      "flow_scale": 1e-3,
      "pressure_scale": 1e5
    }
  }
}
```

Hydraulics works in both modes:
- **`dynamic`**: the network is re-solved at every timestep and tank inventories
  integrate the network-determined flows.
- **`steady`**: the network is solved once; unit logic then evaluates T and
  composition at those pressures/flows. (`pf run` routes steady hydraulic
  flowsheets to the native engine instead of the equation-oriented solver.)

When enabled:
- Stream **pressures** are solved simultaneously for all streams.
- Element **flows** are determined by the pressure–flow laws of valves, pipes,
  strainers, and pumps.
- **Tank** `outlet_flow` is **ignored**: the network determines the actual outflow.
  A nonzero net tank flow is accumulation/depletion in dynamic mode.
- **Pipes** and **valves** must specify hydraulic parameters (`Cv`, `valve_opening`,
  `design_flowrate`, etc.) instead of relying solely on `delta_p`.
- **Pumps** can use `pump_curve` for realistic head-flow behavior or `delta_p` as a
  fixed-rise fallback.
- Declared feed streams that feed an element directly are network boundaries, so
  their flow is solved. Feeds that feed a `Tank` keep their declared `flowrate`
  (the tank inventory integrates it).

### Solver overview

The network solver builds a nonlinear system:
- **Nodes** are streams connected to tanks or declared as feeds are fixed-pressure
  boundaries; every other stream is a free node with a mass balance.
- **Terminal outlets** — a produced stream that nothing consumes, such as a drain
  or product — are anchored at `outlet_pressure` on the producing unit, a declared
  stream pressure, or atmospheric pressure (`101325 Pa`) as a last resort.
- **Edges** are units (Pump, Valve, Pipes, Strainer, Heater). Each edge contributes
  one pressure–flow relation.
- **Junctions** — a unit with multiple inlets (for example a mixing valve) — merge
  its inlet streams into one zero-volume node. The node has one pressure unknown
  and its mass balance sums every incoming branch. A `Pump` must have exactly one
  inlet.

The system is solved with `scipy.optimize.least_squares` using pressure (Pa)
and flow (mol/s) scaling. The flow scale is inferred from the design flow and feed
flowrates unless `flow_scale` is set. An active-set loop closes check valves and
pumps that would reverse flow. On failure, `HydraulicSolveError` is raised with the
residual norm, the worst equation, the solver status, and the active valve/pump
states.

### Property resolution

Hydraulic laws need `density`, `viscosity`, `molar_volume`, and `cp`. These are
resolved from the attached thermodynamic provider (CoolProp by default) or from
explicit unit parameters (`density`, `viscosity`, `molar_volume`, `cp` in the
unit config). Explicit parameters always override provider values.

### Flash in hydraulic mode

`Flash` is **not supported** in hydraulic mode because it has multiple outlets
and would need internal composition-dependent pressure splits. Use hydraulic
mode for networks built from pumps, valves, pipes, strainers, heaters and tanks;
multi-inlet junction units are supported, but multi-outlet units are not.

## Examples

### Default SM units (no provider declared)

```json
{
  "units": {
    "pump_1": { "type": "Pump", "in": "feed", "out": "after_pump",
                "delta_p": 200000, "efficiency": 0.75, "material": 3 },
    "valve_1": { "type": "Valve", "in": "after_pump", "out": "product",
                 "pressure_ratio": 0.5, "material": 3 }
  }
}
```

### Hydraulic mode: closed-loop water system

```json
{
  "simulation": {
    "mode": "dynamic",
    "t0": 0.0, "tf": 20.0, "dt": 1.0,
    "hydraulics": { "enabled": true, "tol": 1e-8, "max_iter": 200 }
  },
  "units": {
    "tank_1": {
      "type": "Tank", "in": ["after_pipe_5"], "out": ["tank_1_to_pipe1"],
      "P": 101325, "initial_n": { "Water": 222000 }, "initial_T": 298.15,
      "material": 1
    },
    "pump_1": {
      "type": "Pump", "in": "after_pipe_1", "out": "after_pump",
      "pump_curve": {
        "flow": [0.0, 0.0005556, 0.0011112],
        "delta_p": [4000000.0, 3000000.0, 1500000.0]
      },
      "material": 2
    },
    "makeup_valve": {
      "type": "Valve", "in": "after_pipe_2", "out": "after_valve",
      "Cv": 1.04e-6, "valve_opening": 0.5, "valve_function": "linear",
      "material": 2
    },
    "return_pipe": {
      "type": "Pipes", "in": "after_valve", "out": "after_pipe_3",
      "delta_p": 1000.0, "diameter": 0.05, "length": 1.0,
      "design_flowrate": 0.0005556, "material": 2
    },
    "tank_2": {
      "type": "Tank", "in": ["after_pipe_3"], "out": "after_tank_2",
      "P": 1961330, "initial_n": { "Water": 555000 }, "initial_T": 298.15,
      "material": 1
    }
  }
}
```

### A unit specifying to the IDAES provider

```json
{
  "units": {
    "pump_1": { "type": "Pump", "provider": "idaes", "in": "feed",
                "out": "after_pump", "delta_p": 200000, "efficiency": 0.75,
                "material": 3 }
  }
}
```
