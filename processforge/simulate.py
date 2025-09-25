import json
import csv
import matplotlib.pyplot as plt
import argparse
from processforge.flowsheet import Flowsheet

def save_results_csv(results, fname="results.csv"):
    """Save steady-state results to a CSV file."""
    comps = set()
    for s in results.values():
        comps.update(s.get("z", {}).keys())
    comps = sorted(comps)

    with open(fname, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["Stream", "T [K]", "P [Pa]", "Phase", "VaporFrac"] + comps
        writer.writerow(header)
        for name, s in results.items():
            row = [
                name,
                s.get("T", ""),
                s.get("P", ""),
                s.get("phase", ""),
                s.get("beta", "")
            ]
            for c in comps:
                row.append(s.get("z", {}).get(c, 0.0))
            writer.writerow(row)

def save_timeseries_csv(results, fname="results_timeseries.csv"):
    """
    Save dynamic results (time series).
    results: dict {time: {stream: {...}}}
    """
    # Collect unique components
    times = sorted(results.keys())
    streams = set()
    comps = set()
    for t in times:
        for sname, sdata in results[t].items():
            streams.add(sname)
            comps.update(sdata.get("z", {}).keys())
    streams = sorted(streams)
    comps = sorted(comps)

    with open(fname, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["time", "stream", "T [K]", "P [Pa]", "Phase", "VaporFrac"] + comps
        writer.writerow(header)
        for t in times:
            for sname in streams:
                s = results[t].get(sname, {})
                row = [
                    t,
                    sname,
                    s.get("T", ""),
                    s.get("P", ""),
                    s.get("phase", ""),
                    s.get("beta", "")
                ]
                for c in comps:
                    row.append(s.get("z", {}).get(c, 0.0))
                writer.writerow(row)

def plot_results(results, fname="results.png"):
    """Generate steady-state charts."""
    streams = list(results.keys())
    temps = [results[s].get("T", None) for s in streams]

    # Temperature bar chart
    plt.figure(figsize=(8, 4))
    plt.bar(streams, temps, color="skyblue")
    plt.ylabel("Temperature [K]")
    plt.title("Stream Temperatures (Steady State)")
    plt.tight_layout()
    plt.savefig("temps_" + fname)
    plt.close()

    # Composition stacked bar chart
    comps = set()
    for s in results.values():
        comps.update(s.get("z", {}).keys())
    comps = sorted(comps)

    bottom = [0.0] * len(streams)
    plt.figure(figsize=(8, 4))
    for comp in comps:
        vals = [results[s].get("z", {}).get(comp, 0.0) for s in streams]
        plt.bar(streams, vals, bottom=bottom, label=comp)
        bottom = [bottom[i] + vals[i] for i in range(len(vals))]
    plt.ylabel("Mole Fraction")
    plt.title("Stream Compositions (Steady State)")
    plt.legend()
    plt.tight_layout()
    plt.savefig("comps_" + fname)
    plt.close()

def plot_timeseries(results, fname="timeseries.png"):
    """Generate time-series plots for dynamic simulations."""
    times = sorted(results.keys())
    streams = list(next(iter(results.values())).keys())

    # Temperature vs time
    plt.figure(figsize=(8, 5))
    for sname in streams:
        T_series = [results[t][sname].get("T", None) for t in times]
        plt.plot(times, T_series, label=sname)
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [K]")
    plt.title("Stream Temperatures vs Time")
    plt.legend()
    plt.tight_layout()
    plt.savefig("temps_" + fname)
    plt.close()

    # Composition vs time (line chart for each component of each stream)
    comps = set()
    for t in times:
        for s in results[t].values():
            comps.update(s.get("z", {}).keys())
    comps = sorted(comps)

    for sname in streams:
        plt.figure(figsize=(8, 5))
        for comp in comps:
            comp_series = [results[t][sname].get("z", {}).get(comp, 0.0) for t in times]
            plt.plot(times, comp_series, label=comp)
        plt.xlabel("Time [s]")
        plt.ylabel("Mole Fraction")
        plt.title(f"Compositions vs Time ({sname})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"comps_{sname}_" + fname)
        plt.close()

def main():
    parser = argparse.ArgumentParser(
        description="Process Forge - Chemical Process Simulation",
        prog="processforge"
    )
    parser.add_argument(
        "flowsheet",
        help="Path to the flowsheet JSON file"
    )

    args = parser.parse_args()

    fname = args.flowsheet
    with open(fname) as f:
        data = json.load(f)
    fs = Flowsheet(data)
    results = fs.run()

    # Detect steady vs dynamic (placeholder: if dict-of-dicts keyed by time, it's dynamic)
    if isinstance(next(iter(results.values())), dict) and "T" in next(iter(results.values())).keys():
        # Steady state
        print("=== Steady-State Results ===")
        for sname, sdata in results.items():
            print(sname, sdata)
        save_results_csv(results, "results.csv")
        plot_results(results, "results.png")
        print("Saved results.csv, temps_results.png, comps_results.png")
    else:
        # Dynamic
        print("=== Dynamic Results ===")
        for t, streams in results.items():
            print(f"t={t}: {streams}")
        save_timeseries_csv(results, "results_timeseries.csv")
        plot_timeseries(results, "timeseries.png")
        print("Saved results_timeseries.csv + time series plots")

if __name__ == "__main__":
    main()
