# Unit Operations Reference

This page documents every unit type that can appear in a processforge flowsheet's `units` block, with the configurable **properties** — and which properties are computed by the **built-in (SM)** logic versus the **IDAES provider** when a unit declares `"provider": "idaes"`.

## 

Processforge deliberately keeps **one unit type name with two property vectors**:

- **Built-in (SM)** — the default sequential-modular models in
  `processforge.units`. 
- **IDAES (provider)** — when a unit sets `"provider": "idaes"`, the IDAES
  provider intercepts the same unit types and can compute some of them with IDAES steady-state models.

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
mole-fraction dict (mutually exclusive with `material_mix`). The built-in
models assume water-like properties (ρ = 1000 kg/m³, Cp = 4180 J/kg·K,
MW = 0.018 kg/mol) unless a thermodynamic provider is attached.

## Unit operation properties

The **Built-in (SM)** column says whether the default logic reads the key; the **IDAES (provider)** column says whether the IDAES provider understands it. A `—`
means the vector does not model that property. Cells marked **fall-through**
are unit types IDAES does **not** intercept — those always run the built-in
logic today, but are reserved for future IDAES features.

Legend:
- ✅ — property is read on this vector.
- — — property is not modelled on this vector.
- 🔽 — unit type falls through to the built-in logic (IDAES does not intercept it).

### Pump

Adds a fixed pressure rise `delta_p`; computes shaft power and outlet
temperature.

| Property | Built-in (SM) | IDAES (provider) | Default | Units | Description |
|----------|---------------|------------------|---------|-------|-------------|
| `delta_p` | ✅ | ✅ | `1e5` | Pa | Pressure rise. |
| `efficiency` | ✅ | ✅ | `0.8` | — | Adiabatic efficiency. |

### Valve

Isenthalpic pressure reduction by a fixed ratio.

| Property | Built-in (SM) | IDAES (provider) | Default | Units | Description |
|----------|---------------|------------------|---------|-------|-------------|
| `pressure_ratio` | ✅ | ✅ | `0.5` | — | Outlet pressure / inlet pressure. |

### Strainer

Fixed pressure drop element.

| Property | Built-in (SM) | IDAES (provider) | Default | Units | Description |
|----------|---------------|------------------|---------|-------|-------------|
| `delta_p` | ✅ | ✅ | `5000.0` | Pa | Pressure drop. |

### Pipes

Laminar pipe with a fixed pressure drop; optionally computes flow from the
drop via a simplified Hagen–Poiseuille relation.

| Property | Built-in (SM) | IDAES (provider) | Default | Units | Description |
|----------|---------------|------------------|---------|-------|-------------|
| `delta_p` | ✅ | ✅ | `1000.0` | Pa | Pressure drop. |
| `diameter` | ✅ | ✅ | `0.1` | m | Pipe diameter. |

### Tank

Well-mixed molar tank. **Dynamic only** (requires dynamic mode).

| Property | Built-in (SM) | IDAES (provider) | Default | Units | Description |
|----------|---------------|------------------|---------|-------|-------------|
| `outlet_flow` | ✅ | — | `1.0` | mol/s | Outlet molar flow. |
| `P` | ✅ | — | `101325.0` | Pa | Operating pressure. |
| `duty` | ✅ | — | `0.0` | W | Heat duty (positive heats). |
| `initial_n` | ✅ | — | `{}` | mol | Initial molar inventory. |
| `initial_T` | ✅ | — | `300.0` | K | Initial temperature. |
| `flow_in` | ✅ | — | `1.0` | mol/s | Inlet flow when stream lacks one. |

### Flash

Isothermal flash separator (Rachford–Rice) producing vapor and liquid
outlets. 🔽 **Fall-through** — needs K-values from a thermo provider.**

| Property | Built-in (SM) | IDAES (provider) | Default | Units | Description |
|----------|---------------|------------------|---------|-------|-------------|
| `P` | ✅ | — | — | Pa | Operating pressure. |
| `in` | ✅ | — | — | — | Inlet stream name. |
| `out_vap` | ✅ | — | — | — | Outlet vapor stream name. |
| `out_liq` | ✅ | — | — | — | Outlet liquid stream name. |

## Heater, CSTR and PFR

These units are not supported in IDAES currently. 


## Provider resolution & fall-through

IDAES does **not** define new unit classes. It intercepts the same SM unit
names via `IdaesProvider.compute_unit()`:

| Unit type | IDAES intercepts? | IDAES reads | Notes |
|-----------|-------------------|-------------|-------|
| `Pump` | ✅ | `delta_p`, `efficiency` | — |
| `Valve` | ✅ | `pressure_ratio` | — |
| `Strainer` | ✅ | `delta_p` | — |
| `Pipes` | ✅ | `delta_p`, `diameter` | — |
| `Heater` | 🔽 | — | Needs thermo provider; built-in logic runs. |
| `Flash` | 🔽 | — | Needs K-values; built-in logic runs. |
| `CSTR` / `PFR` / `Tank` | 🔽 | — | Default SM / Cantera logic runs. |

Execution order per unit: built-in SM logic → provider (`compute_unit`) →
fall-through to the built-in SM logic when the provider returns `None`.

## Examples

Default SM units (no provider declared):

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

A unit specifying to the IDAES provider:

```json
{
  "units": {
    "pump_1": { "type": "Pump", "provider": "idaes", "in": "feed",
                "out": "after_pump", "delta_p": 200000, "efficiency": 0.75,
                "material": 3 }
  }
}
```
