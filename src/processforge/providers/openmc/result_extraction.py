"""OpenMC statepoint result extraction helpers."""
from __future__ import annotations

import math
import pathlib

from loguru import logger

from processforge.schemas.openmc.openmc_model import OpenMCSetting, RunMode
from processforge.types import OutputArtifact, OutputField, Quantity


# Score → SI-ish unit string, so outputs are comparable with other engines
# (CoolProp J/mol, FESTIM W/m², …) via the shared ``Quantity``/pint layer.
_SCORE_UNITS = {
    "flux": "n/cm^2/s",
    "fission": "1/cm^3/s",
    "heating": "W",
    "kappa_fission": "W",
    "absorption": "1/cm^3/s",
    "scatter": "1/cm^3/s",
    "elastic": "1/cm^3/s",
    "capture": "1/cm^3/s",
    "total": "1/cm^3/s",
    "current": "1/cm^2/s",
    "delayed_nu_fission": "1/cm^3/s",
}
# Energy per fission used to convert a fission rate into a comparable power.
_MEV_PER_FISSION = 200.0
_J_PER_MEV = 1.602176634e-13  # 1 MeV = 1e6 eV * 1.602e-19 J/eV


def extract_results(
    openmc,
    solver_cfg: OpenMCSetting,
    tally_cfgs: list,
    statepoint_path,
    run_dir: pathlib.Path,
) -> tuple:
    """Parse the statepoint and return ``(fields, artifacts, diagnostics)``.

    Standardized, unit-bearing fields
    ----------------------------------
    * ``k_eff`` (dimensionless) with ``std_dev`` — eigenvalue mode only
    * Per mesh tally + score:
      ``tally_{id}_{score}_mean``       — volume-averaged bin mean
      ``tally_{id}_{score}_integrated`` — Σ(bin mean × cell volume)
      both tagged with the score's unit (e.g. ``flux`` → ``n/cm^2/s``).
    * ``power`` (W) — derived from the total fission rate
      (``fission`` score) using 200 MeV/fission, so OpenMC outputs are
      directly comparable with CoolProp/FESTIM thermal quantities.

    The raw per-voxel tally dataframe is also written to a CSV
    :class:`OutputArtifact`, and the statepoint HDF5 is registered as an
    artifact. Extraction issues are collected in ``diagnostics["tally_warnings"]``.
    """
    fields: list = []
    artifacts: list = []
    diagnostics: dict = {"run_dir": str(run_dir.resolve())}
    warnings: list = []

    if statepoint_path is None:
        warnings.append("run produced no statepoint file")
        diagnostics["tally_warnings"] = warnings
        return fields, artifacts, diagnostics

    sp_path = str(pathlib.Path(statepoint_path).resolve())
    artifacts.append(OutputArtifact(
        name="statepoint",
        kind="statepoint",
        local_path=sp_path,
        source="local",
    ))
    diagnostics["statepoint_path"] = sp_path

    sp = openmc.StatePoint(statepoint_path)
    try:
        # k_eff — only present in eigenvalue mode
        if sp.keff is not None:
            fields.append(OutputField(
                name="k_eff",
                quantity=Quantity(
                    value=float(sp.keff.n), unit="", std_dev=float(sp.keff.s)
                ),
                kind="scalar",
                source="keff",
            ))
            logger.info(
                f"OpenMCProvider: k_eff = {float(sp.keff.n):.6f} "
                f"+/- {float(sp.keff.s):.6f}"
            )
        elif solver_cfg.run_mode == RunMode.eigenvalue:
            warnings.append("eigenvalue run produced no k_eff in statepoint")

        for tally_cfg in tally_cfgs:
            try:
                tally = sp.get_tally(id=tally_cfg.tally_id)
            except Exception as exc:  # noqa: BLE001
                warnings.append(
                    f"tally id={tally_cfg.tally_id} not found in statepoint: {exc}"
                )
                continue

            # Cell volume (area for 2D meshes) for volume-weighted aggregation.
            try:
                ll = [float(x) for x in tally.mesh.lower_left]
                ur = [float(x) for x in tally.mesh.upper_right]
                dim = [float(x) for x in tally.mesh.dimension]
                cell_vol = 1.0
                for a, b, d in zip(ll, ur, dim):
                    cell_vol *= (b - a) / d if d else 1.0
            except Exception:  # noqa: BLE001
                cell_vol = 1.0

            for score in tally_cfg.scores:
                try:
                    df = tally.get_pandas_dataframe(scores=[score])
                    means = [float(x) for x in df["mean"].values]
                    std_devs = [float(x) for x in df["std. dev."].values]
                    unit = _SCORE_UNITS.get(score, "")
                    key_prefix = f"tally_{tally_cfg.tally_id}_{score}"

                    mean_val = sum(means) / len(means) if means else 0.0
                    integrated = sum(m * cell_vol for m in means)
                    integrated_std = math.sqrt(
                        sum((s * cell_vol) ** 2 for s in std_devs)
                    ) if std_devs else 0.0

                    fields.append(OutputField(
                        name=f"{key_prefix}_mean",
                        quantity=Quantity(value=mean_val, unit=unit),
                        kind="scalar",
                        source=f"tally_{tally_cfg.tally_id}/{score}",
                    ))
                    fields.append(OutputField(
                        name=f"{key_prefix}_integrated",
                        quantity=Quantity(
                            value=integrated, unit=unit, std_dev=integrated_std
                        ),
                        kind="scalar",
                        source=f"tally_{tally_cfg.tally_id}/{score}",
                    ))

                    # Total fission rate → comparable power (W).
                    if score == "fission" and integrated > 0:
                        power_W = (
                            integrated
                            * _MEV_PER_FISSION
                            * _J_PER_MEV
                        )
                        power_std = (
                            power_W * (integrated_std / integrated)
                            if integrated_std else 0.0
                        )
                        fields.append(OutputField(
                            name="power",
                            quantity=Quantity(
                                value=power_W, unit="W", std_dev=power_std
                            ),
                            kind="scalar",
                            source="fission_rate",
                        ))
                        diagnostics.setdefault("notes", []).append(
                            "power derived from total fission rate assuming "
                            f"{_MEV_PER_FISSION} MeV/fission."
                        )
                except Exception as exc:  # noqa: BLE001
                    warnings.append(
                        f"could not extract score '{score}' from tally "
                        f"{tally_cfg.tally_id}: {exc}"
                    )

            # Raw per-voxel field as a CSV artifact for downstream post-processing.
            try:
                csv_name = f"tally_{tally_cfg.tally_id}_{tally_cfg.name or tally_cfg.tally_id}.csv"
                csv_path = pathlib.Path(run_dir) / csv_name
                df_all = tally.get_pandas_dataframe()
                df_all.to_csv(csv_path)
                artifacts.append(OutputArtifact(
                    name=csv_name,
                    kind="csv",
                    local_path=str(csv_path.resolve()),
                    source="local",
                ))
            except Exception as exc:  # noqa: BLE001
                warnings.append(
                    f"could not export tally {tally_cfg.tally_id} CSV: {exc}"
                )
    finally:
        del sp

    if warnings:
        diagnostics["tally_warnings"] = warnings
    return fields, artifacts, diagnostics
