"""Shared helpers for flowsheet topology inspection."""

TOPOLOGY_KEYS = {"type", "in", "out", "out_liq", "out_vap"}


def get_inlets(unit_cfg: dict) -> list[str]:
    """Return inlet stream name(s) as a list."""
    inlet = unit_cfg.get("in")
    if isinstance(inlet, list):
        return inlet
    return [inlet] if inlet else []


def get_outlets(unit_cfg: dict) -> list[str]:
    """Return outlet stream name(s), including Flash vapor/liquid branches."""
    if unit_cfg.get("type") == "Flash":
        outlets = []
        if unit_cfg.get("out_liq"):
            outlets.append(unit_cfg["out_liq"])
        if unit_cfg.get("out_vap"):
            outlets.append(unit_cfg["out_vap"])
        return outlets
    out = unit_cfg.get("out")
    if isinstance(out, list):
        return out
    return [out] if out else []
