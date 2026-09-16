"""Unified ProcessStateArchive — the single store for simulation state + outputs.

A ``<base>.pfarchive`` directory holds, in one place:

* ``snapshots/``         — converged EO state vectors (delegated to the proven
  :class:`~processforge.state.StateManager`, preserving drift detection,
  homotopy deltas and rollback).
* ``runs/<run_id>.json`` — a :class:`~processforge.types.RunManifest` capturing
  every engine output (OpenMC, FESTIM, CoolProp, Cantera) for one ``pf run`` /
  ``pf apply``.
* ``results/<run_id>/results.zarr`` — stream timeseries and engine outputs
  in a compressed Zarr store (the single source of truth for time-series data).
* ``latest_run``         — pointer to the most recent run.
"""
from __future__ import annotations

import json
import os
from typing import Optional

from loguru import logger

from ..state import StateManager
from ..types import RunManifest, SnapshotState


def _read_json(path: str, default):
    if not os.path.exists(path):
        return default
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def _write_json(path: str, data) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


class ProcessStateArchive:
    """Single archive for a flowsheet's state + standardized engine outputs."""

    def __init__(self, path: str) -> None:
        self.path = path
        os.makedirs(self.path, exist_ok=True)
        self._state = StateManager(self.path)

    # ==================================================================
    # Snapshot (state-vector) delegation
    # ==================================================================
    def save_snapshot(self, config, x_converged, var_names, *, metadata=None, parent_snapshot_id=None) -> str:
        return self._state.save_state(
            config, x_converged, var_names,
            metadata=metadata, parent_snapshot_id=parent_snapshot_id,
        )

    def load_snapshot(self) -> Optional[SnapshotState]:
        return self._state.load_state()

    def detect_drift(self, current_config, state) -> list[str]:
        return self._state.detect_drift(current_config, state)

    def detect_structural_diff(self, current_config, state) -> dict:
        return self._state.detect_structural_diff(current_config, state)

    def validate_metadata(self, current_metadata, state) -> list[str]:
        return self._state.validate_metadata(current_metadata, state)

    def rollback(self, n: int = 1) -> bool:
        return self._state.rollback(n)

    def state_to_stream_dicts(self, state) -> dict:
        return self._state.state_to_stream_dicts(state)

    def list_snapshots(self) -> list[dict]:
        return self._state.list_snapshots()

    # ==================================================================
    # Runs (engine outputs)
    # ==================================================================
    def _runs_dir(self) -> str:
        return os.path.join(self.path, "runs")

    def save_run(self, manifest: RunManifest) -> str:
        """Persist a :class:`RunManifest` to disk."""
        os.makedirs(self._runs_dir(), exist_ok=True)
        run_path = os.path.join(self._runs_dir(), f"{manifest.run_id}.json")
        _write_json(run_path, json.loads(manifest.model_dump_json()))

        with open(os.path.join(self.path, "latest_run"), "w", encoding="utf-8") as f:
            f.write(manifest.run_id)

        logger.info(f"Saved run '{manifest.run_id}' to {self.path}")
        return manifest.run_id

    def load_run(self, run_id: str) -> Optional[RunManifest]:
        run_path = os.path.join(self._runs_dir(), f"{run_id}.json")
        if not os.path.exists(run_path):
            return None
        with open(run_path, "r", encoding="utf-8") as f:
            return RunManifest(**json.load(f))

    def latest_run(self) -> Optional[RunManifest]:
        latest = os.path.join(self.path, "latest_run")
        if not os.path.exists(latest):
            return None
        with open(latest, encoding="utf-8") as f:
            run_id = f.read().strip()
        return self.load_run(run_id)

    def list_runs(self) -> list[str]:
        if not os.path.isdir(self._runs_dir()):
            return []
        return sorted(
            p.replace(".json", "")
            for p in os.listdir(self._runs_dir())
            if p.endswith(".json")
        )

    def field_occurrences(self, field_name: str) -> list[dict]:
        """Scan ``runs/*.json`` for fields matching *field_name*."""
        results: list[dict] = []
        runs_dir = self._runs_dir()
        if not os.path.isdir(runs_dir):
            return results
        for fname in sorted(os.listdir(runs_dir)):
            if not fname.endswith(".json"):
                continue
            run_id = fname.replace(".json", "")
            manifest = self.load_run(run_id)
            if manifest is None:
                continue
            for unit_name, out in manifest.units.items():
                for f in out.fields:
                    if f.name == field_name:
                        results.append({
                            "run_id": run_id, "scope": "unit",
                            "name": unit_name, "engine": out.engine,
                            "units": f.quantity.unit,
                        })
            for stream_name, out in manifest.streams.items():
                for f in out.fields:
                    if f.name == field_name:
                        results.append({
                            "run_id": run_id, "scope": "stream",
                            "name": stream_name, "engine": out.engine,
                            "units": f.quantity.unit,
                        })
        return results
