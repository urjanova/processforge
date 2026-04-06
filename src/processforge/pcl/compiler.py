"""PCL compiler: load a .pcl (TOML) file and transform it to a flowsheet config dict."""
import tomllib
from .namespace import PCL_NAMESPACE


class PCLCompileError(ValueError):
    pass


def load_pcl(path: str) -> dict:
    """Load a .pcl file and return a compiled flowsheet config dict."""
    with open(path, "rb") as f:
        raw = tomllib.load(f)
    return compile_to_dict(raw, source=path)


def compile_to_dict(raw: dict, source: str = "<pcl>") -> dict:
    """
    Post-transform a tomllib-loaded PCL dict into a flowsheet config dict.

    Transforms applied:
    - units[name].source = "pf.*.*"  →  units[name].type = "<UnitType>" (source key removed)
    - streams[name].units_annotation  →  streams[name]._units (key renamed)
    """
    result = dict(raw)

    # Transform units
    if "units" in result:
        transformed_units = {}
        for unit_name, unit_cfg in result["units"].items():
            unit = dict(unit_cfg)
            if "source" in unit:
                ref = unit.pop("source")
                if ref not in PCL_NAMESPACE:
                    raise PCLCompileError(
                        f"Unknown provider reference '{ref}' at units.{unit_name!r} "
                        f"in {source!r}. "
                        f"Valid references: {sorted(PCL_NAMESPACE)}"
                    )
                unit["type"] = PCL_NAMESPACE[ref]
            transformed_units[unit_name] = unit
        result["units"] = transformed_units

    # Transform stream units_annotation → _units
    if "streams" in result:
        transformed_streams = {}
        for stream_name, stream_cfg in result["streams"].items():
            stream = dict(stream_cfg)
            if "units_annotation" in stream:
                stream["_units"] = stream.pop("units_annotation")
            transformed_streams[stream_name] = stream
        result["streams"] = transformed_streams

    return result
