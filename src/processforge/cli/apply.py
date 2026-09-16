"""``pf apply`` — state-based warm start with homotopy fallback."""

from __future__ import annotations

from typing import Literal

import typer
from loguru import logger

from ..runner import ProcessforgeRunError, apply_flowsheet


def apply(
    flowsheet: str = typer.Argument(help="Path to the flowsheet JSON file"),
    backend: Literal["scipy", "pyomo", "casadi"] | None = typer.Option(
        None,
        "--backend",
        help="Override the flowsheet's simulation.backend",
    ),
    tolerance: float = typer.Option(
        1e-6,
        "--tolerance",
        help="Newton solver convergence tolerance (default: 1e-6)",
    ),
    max_iter: int = typer.Option(
        50,
        "--max-iter",
        help="Max Newton iterations (default: 50)",
    ),
    skip_homotopy: bool = typer.Option(
        False,
        "--skip-homotopy",
        help="Disable homotopy fallback; cold-start only",
    ),
) -> None:
    """Apply flowsheet: drift detection, warm-start, homotopy fallback, convergence guardrails."""
    try:
        result = apply_flowsheet(
            flowsheet,
            backend=backend,
            tolerance=tolerance,
            max_iter=max_iter,
            skip_homotopy=skip_homotopy,
        )
    except ProcessforgeRunError as exc:
        logger.error(str(exc))
        raise SystemExit(1) from exc

    logger.info("=== Apply Summary ===")
    logger.info(f"  Status      : {result.status}")
    logger.info(f"  Run ID      : {result.run_id}")
    logger.info(f"  Snapshot ID : {result.snapshot_id}")
    logger.info(f"  Archive     : {result.archive_path}")
    if result.remote_uris:
        logger.info(f"  Remote URIs : {len(result.remote_uris)} object(s) in S3")
