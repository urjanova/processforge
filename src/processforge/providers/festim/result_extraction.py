"""FESTIM result extraction helpers."""
from __future__ import annotations

import csv
import pathlib

import numpy as np
from loguru import logger

from processforge.types import Axis, Domain, OutputArtifact, OutputField, Quantity


def extract_results(exports, run_dir: pathlib.Path) -> tuple:
    """Return ``(fields, artifacts, diagnostics)`` from FESTIM exports.

    * Every export exposing a ``title`` + ``value`` (surface fluxes, total
      volumes) contributes a scalar :class:`OutputField`.
    * Every export exposing ``t`` + ``data`` (incl. ``Profile1DExport``)
      contributes a reduced ``{title}_mean_total`` / ``{title}_std_dev``
      :class:`OutputField` (mirroring the OpenMC tally convention) plus the
      full series written to a content-addressed CSV :class:`OutputArtifact`.
      The field carries a ``timeseries`` ``Domain`` so downstream code can
      reconstruct the shape without the file.

    ``diagnostics`` always carries ``run_dir``.
    """
    fields: list = []
    artifacts: list = []
    diagnostics: dict = {"run_dir": str(run_dir.resolve())}

    for export in exports:
        title = getattr(export, "title", None)
        value = getattr(export, "value", None)
        if title is not None and value is not None:
            fields.append(OutputField(
                name=str(title),
                quantity=Quantity(value=[float(value)], unit=""),
                kind="scalar",
                source=str(title),
            ))

        t = getattr(export, "t", None)
        data = getattr(export, "data", None)
        if not (t and data):
            continue

        key = str(title if title is not None else getattr(export, "field", ""))
        if not key:
            continue

        # Flatten all time steps into a single value array for statistics.
        flat: list[float] = []
        rows: list[list[float]] = []
        for step in data:
            if hasattr(step, "tolist"):
                arr = step.tolist()
            else:
                arr = [float(step)]
            if not isinstance(arr, list):
                arr = [arr]
            rows.append(arr)
            flat.extend(arr)

        flat_arr = np.asarray(flat, dtype=float)
        mean_total = float(flat_arr.mean()) if flat_arr.size else 0.0
        std_dev = float(flat_arr.std()) if flat_arr.size else 0.0

        domain = Domain(
            type="timeseries",
            axes=[Axis(name="t", unit="s", coordinates=list(map(float, t)))],
            shape=[len(rows), max((len(r) for r in rows), default=0)],
        )
        fields.append(OutputField(
            name=f"{key}_mean_total",
            quantity=Quantity(value=mean_total, unit=""),
            kind="timeseries",
            source=f"{key}/mean",
            domain=domain,
        ))
        fields.append(OutputField(
            name=f"{key}_std_dev",
            quantity=Quantity(value=std_dev, unit=""),
            kind="timeseries",
            source=f"{key}/std_dev",
            domain=domain,
        ))

        # Persist the full series as a CSV artifact (not inlined over HTTP).
        csv_name = f"{key}.csv"
        csv_path = pathlib.Path(run_dir) / csv_name
        try:
            with open(csv_path, "w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh)
                writer.writerow(["t"] + [f"v{i}" for i in range(max((len(r) for r in rows), default=0))])
                for ti, row in zip(t, rows):
                    writer.writerow([float(ti)] + list(row))
            artifacts.append(OutputArtifact(
                name=csv_name,
                kind="csv",
                local_path=str(csv_path.resolve()),
                source="local",
            ))
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"FestimProvider: failed to write series CSV '{csv_name}': {exc}")

    return fields, artifacts, diagnostics
