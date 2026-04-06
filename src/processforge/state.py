"""StateManager: versioned .pfstate Zarr snapshot store with drift detection."""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone

import numpy as np
import zarr
from loguru import logger


class StateManager:
    """
    Manages the ``.pfstate`` versioned snapshot store.

    Layout::

        outputs/my-flowsheet.pfstate/
            snapshots/
                0001_2026-04-06T12:00:00Z/   ← Zarr group per successful apply
                    .attrs  → {config, var_names, timestamp, snapshot_id}
                    x       ← float64 dataset
                0002_2026-04-06T13:30:00Z/
                    ...
            latest                            ← plain-text pointer: snapshot dir name

    Args:
        state_path: Path to the ``.pfstate`` directory (e.g. ``outputs/my.pfstate``).
    """

    def __init__(self, state_path: str) -> None:
        self.state_path = state_path
        self._snapshots_dir = os.path.join(state_path, "snapshots")
        self._latest_file = os.path.join(state_path, "latest")

    # ------------------------------------------------------------------
    # Snapshot persistence
    # ------------------------------------------------------------------

    def save_state(
        self,
        config: dict,
        x_converged: np.ndarray,
        var_names: list[str],
    ) -> str:
        """Create a new numbered snapshot and update the ``latest`` pointer.

        Returns:
            The snapshot directory name (e.g. ``"0002_2026-04-06T13:30:00Z"``).
        """
        os.makedirs(self._snapshots_dir, exist_ok=True)

        # Determine next snapshot number
        existing = self._list_snapshot_dirs()
        next_num = len(existing) + 1
        ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        snapshot_name = f"{next_num:04d}_{ts}"
        snapshot_path = os.path.join(self._snapshots_dir, snapshot_name)

        store = zarr.storage.LocalStore(snapshot_path)
        root = zarr.group(store=store, overwrite=True)
        root.attrs["config"] = json.dumps(config)
        root.attrs["var_names"] = list(var_names)
        root.attrs["timestamp"] = ts
        root.attrs["snapshot_id"] = snapshot_name
        
        x_data = np.asarray(x_converged, dtype=float)
        root.create_dataset("x", data=x_data, shape=x_data.shape)

        # Update latest pointer
        with open(self._latest_file, "w", encoding="utf-8") as f:
            f.write(snapshot_name)

        logger.info(f"Saved snapshot {snapshot_name} to {self.state_path}")
        return snapshot_name

    def load_state(self) -> dict | None:
        """Load the latest snapshot.  Returns ``None`` if no state exists."""
        if not os.path.exists(self._latest_file):
            return None

        with open(self._latest_file, encoding="utf-8") as f:
            latest_name = f.read().strip()

        snapshot_path = os.path.join(self._snapshots_dir, latest_name)
        if not os.path.exists(snapshot_path):
            logger.warning(f"Latest snapshot '{latest_name}' not found on disk.")
            return None

        try:
            store = zarr.storage.LocalStore(snapshot_path)
            root = zarr.group(store=store, overwrite=False)
            config = json.loads(root.attrs["config"])
            x = root["x"][:]
            var_names = list(root.attrs["var_names"])
            return {
                "config": config,
                "x": x,
                "var_names": var_names,
                "snapshot_id": root.attrs.get("snapshot_id", latest_name),
                "timestamp": root.attrs.get("timestamp", ""),
            }
        except Exception as exc:
            logger.warning(f"Failed to load snapshot '{latest_name}': {exc}")
            return None

    def list_snapshots(self) -> list[dict]:
        """Return metadata for all snapshots (oldest first)."""
        dirs = self._list_snapshot_dirs()
        snapshots = []
        for d in dirs:
            path = os.path.join(self._snapshots_dir, d)
            try:
                store = zarr.storage.LocalStore(path)
                root = zarr.group(store=store, overwrite=False)
                snapshots.append({
                    "id": root.attrs.get("snapshot_id", d),
                    "timestamp": root.attrs.get("timestamp", ""),
                    "snapshot_dir": path,
                })
            except Exception:
                snapshots.append({"id": d, "timestamp": "", "snapshot_dir": path})
        return snapshots

    def rollback(self, n: int = 1) -> bool:
        """Revert the ``latest`` pointer ``n`` snapshots back.

        Returns:
            ``True`` if rollback succeeded, ``False`` if there is no older snapshot.
        """
        dirs = self._list_snapshot_dirs()
        if not dirs:
            logger.warning("No snapshots to roll back to.")
            return False

        current_name = ""
        if os.path.exists(self._latest_file):
            with open(self._latest_file, encoding="utf-8") as f:
                current_name = f.read().strip()

        try:
            current_idx = dirs.index(current_name)
        except ValueError:
            current_idx = len(dirs) - 1  # fallback: treat as last

        target_idx = current_idx - n
        if target_idx < 0:
            logger.warning("Cannot roll back further — already at the oldest snapshot.")
            return False

        target_name = dirs[target_idx]
        with open(self._latest_file, "w", encoding="utf-8") as f:
            f.write(target_name)
        logger.info(f"Rolled back to snapshot {target_name}")
        return True

    # ------------------------------------------------------------------
    # Drift & structural diff
    # ------------------------------------------------------------------

    def detect_drift(self, current_config: dict, state: dict) -> list[str]:
        """Return parameter paths that differ between ``current_config`` and saved state.

        Paths look like ``"streams.feed.T"`` or ``"units.pump_1.parameters.deltaP"``.
        """
        old_config = state["config"]
        drifted: list[str] = []

        for s_name, s_data in current_config.get("streams", {}).items():
            old_s = old_config.get("streams", {}).get(s_name, {})
            for param in ("T", "P", "flowrate"):
                if s_data.get(param) != old_s.get(param):
                    drifted.append(f"streams.{s_name}.{param}")

            cur_z = s_data.get("z", {})
            old_z = old_s.get("z", {})
            for comp in cur_z:
                if cur_z.get(comp) != old_z.get(comp):
                    drifted.append(f"streams.{s_name}.z.{comp}")

        for u_name, u_cfg in current_config.get("units", {}).items():
            old_u = old_config.get("units", {}).get(u_name, {})
            for param, val in u_cfg.get("parameters", {}).items():
                if val != old_u.get("parameters", {}).get(param):
                    drifted.append(f"units.{u_name}.parameters.{param}")

        return drifted

    def detect_structural_diff(self, current_config: dict, state: dict) -> dict:
        """Compare unit topology between ``current_config`` and saved state.

        Returns a dict with keys ``"added"``, ``"removed"``, ``"modified"``.

        - ``added``: ``{unit_name: unit_type}`` — new units not in old state
        - ``removed``: ``{unit_name: unit_type}`` — units present in old state but not current
        - ``modified``: ``{unit_name: {type, changes}}`` — units present in both but with differing config
        """
        old_units: dict = state["config"].get("units", {})
        new_units: dict = current_config.get("units", {})

        old_names = set(old_units)
        new_names = set(new_units)

        added = {n: new_units[n].get("type", "?") for n in new_names - old_names}
        removed = {n: old_units[n].get("type", "?") for n in old_names - new_names}
        modified: dict = {}

        for name in old_names & new_names:
            old_cfg = old_units[name]
            new_cfg = new_units[name]
            changes: list[str] = []

            # Type change
            if old_cfg.get("type") != new_cfg.get("type"):
                changes.append(f"type: {old_cfg.get('type')} → {new_cfg.get('type')}")

            # Connection changes (in/out)
            for key in ("in", "out", "out_vap", "out_liq"):
                ov, nv = old_cfg.get(key), new_cfg.get(key)
                if ov != nv:
                    changes.append(f"{key}: {ov!r} → {nv!r}")

            # Parameter changes (scalar unit params, not nested "parameters" dict)
            all_keys = set(old_cfg) | set(new_cfg)
            skip = {"type", "in", "out", "out_vap", "out_liq", "provider", "parameters"}
            for k in all_keys - skip:
                ov, nv = old_cfg.get(k), new_cfg.get(k)
                if ov != nv:
                    changes.append(f"{k}: {ov} → {nv}")

            # Nested "parameters" dict
            old_params = old_cfg.get("parameters", {})
            new_params = new_cfg.get("parameters", {})
            for pk in set(old_params) | set(new_params):
                ov, nv = old_params.get(pk), new_params.get(pk)
                if ov != nv:
                    changes.append(f"parameters.{pk}: {ov} → {nv}")

            if changes:
                modified[name] = {"type": new_cfg.get("type", "?"), "changes": changes}

        return {"added": added, "removed": removed, "modified": modified}

    # ------------------------------------------------------------------
    # State → stream dict conversion (for dynamic t=0 initialisation)
    # ------------------------------------------------------------------

    @staticmethod
    def state_to_stream_dicts(state: dict) -> dict[str, dict]:
        """Convert a saved ``x`` vector + ``var_names`` back to stream dicts.

        Returns:
            ``{stream_name: {"T": ..., "P": ..., "flowrate": ..., "z": {...}}}``
        """
        x: np.ndarray = np.asarray(state["x"], dtype=float)
        var_names: list[str] = state.get("var_names", [])

        streams: dict[str, dict] = {}
        for i, label in enumerate(var_names):
            if "/" not in label:
                continue
            stream_name, field = label.split("/", 1)
            if stream_name not in streams:
                streams[stream_name] = {"z": {}}
            if field.startswith("z_"):
                comp = field[2:]
                streams[stream_name]["z"][comp] = float(x[i])
            else:
                streams[stream_name][field] = float(x[i])

        return streams

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _list_snapshot_dirs(self) -> list[str]:
        """Return sorted list of snapshot directory names (oldest first)."""
        if not os.path.isdir(self._snapshots_dir):
            return []
        dirs = sorted(
            d for d in os.listdir(self._snapshots_dir)
            if os.path.isdir(os.path.join(self._snapshots_dir, d))
        )
        return dirs
