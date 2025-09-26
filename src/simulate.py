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
    save_results_json,
    save_timeseries_json,
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
    # Detect steady vs dynamic: if keys are numeric (times), it's dynamic
    is_dynamic = all(isinstance(k, (int, float)) for k in results.keys())
    
    base_name = os.path.splitext(os.path.basename(fname))[0]

    if is_dynamic:
        logger.info("=== Dynamic Results ===")
        save_timeseries_json(results, f"{base_name}_timeseries.json")
        save_results_json(results, f"{base_name}_results.json")  # Assuming this saves summary or something
        save_timeseries_csv(results, f"{base_name}_timeseries.csv")
        plot_timeseries(results, f"{base_name}_timeseries.png")
        logger.info(f"Saved {base_name}_timeseries.csv and {base_name}_timeseries.png")
    else:
        logger.info("=== Steady-State Results ===")
        save_results_json(results, f"{base_name}_results.json")
        save_results_csv(results, f"{base_name}_results.csv")
        save_timeseries_json(results, f"{base_name}_timeseries.json")  # If applicable for steady state
        plot_results(results, f"{base_name}_results.png")
        logger.info(f"Saved {base_name}_results.csv and {base_name}_results.png")

if __name__ == "__main__":
    main()
