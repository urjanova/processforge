import json
import csv
import argparse
from src.flowsheet import Flowsheet
import os
from src.result import (
    save_results_csv,
    save_timeseries_csv,
    plot_results,
    plot_timeseries,
)
from loguru import logger


def main():
    parser = argparse.ArgumentParser(
        description="Process Forge - Chemical Process Simulation", prog="processforge"
    )
    parser.add_argument("flowsheet", help="Path to the flowsheet JSON file")

    args = parser.parse_args()

    fname = args.flowsheet
    if not os.path.exists(fname):
        logger.error(f"Error: Flowsheet file '{fname}' not found.")
        exit(1)
    with open(fname) as f:
        data = json.load(f)
    fs = Flowsheet(data)
    results = fs.run()

    # Detect steady vs dynamic (placeholder: if dict-of-dicts keyed by time, it's dynamic)
    if (
        isinstance(next(iter(results.values())), dict)
        and "T" in next(iter(results.values())).keys()
    ):
        # Steady state
        logger.info("=== Steady-State Results ===")

        save_results_csv(results, "results.csv")
        plot_results(results, "results.png")
        logger.info("Saved results.csv, temps_results.png, comps_results.png")
    else:
        # Dynamic
        logger.info("=== Dynamic Results ===")

        save_timeseries_csv(results, "results_timeseries.csv")
        plot_timeseries(results, "timeseries.png")
        logger.info("Saved results_timeseries.csv + time series plots")


if __name__ == "__main__":
    main()
