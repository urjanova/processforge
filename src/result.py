import csv
import json
import os

import matplotlib.pyplot as plt

def save_results_csv(results, fname="results.csv"):
    """Save steady-state results to a CSV file."""
    comps = set()
    for s in results.values():
        comps.update(s.get("z", {}).keys())
    comps = sorted(comps)

    os.makedirs("outputs", exist_ok=True)
    with open(os.path.join("outputs", fname), "w", encoding="utf-8",newline="") as f:
        writer = csv.writer(f)
        header = ["Stream", "T [K]", "P [Pa]", "Phase", "VaporFrac"] + comps
        writer.writerow(header)
        for name, s in results.items():
            row = [
                name,
                s.get("T", ""),
                s.get("P", ""),
                s.get("phase", ""),
                s.get("beta", ""),
            ]
            for c in comps:
                row.append(s.get("z", {}).get(c, 0.0))
            writer.writerow(row)




def save_timeseries_csv(results, fname="results_timeseries.csv"):
    """
    Save dynamic results (time series).
    results: dict {stream: {prop: [values]}}
    """
    # Collect all streams, components, and the time array
    streams = sorted(results.keys())
    if not streams:
        return
    
    # Assume all streams have the same time array
    times = results[streams[0]].get("time", [])
    
    comps = set()
    for s_data in results.values():
        # The 'z' value can be a dict of lists
        if 'z' in s_data and isinstance(s_data['z'], dict):
            comps.update(s_data['z'].keys())
    comps = sorted(list(comps))

    os.makedirs("outputs", exist_ok=True)
    with open(os.path.join("outputs", fname), "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        header = ["time", "stream", "T [K]", "P [Pa]", "flowrate"] + comps
        writer.writerow(header)

        for i, t in enumerate(times):
            for s_name in streams:
                s_data = results[s_name]
                row = [
                    t,
                    s_name,
                    s_data.get("T", [])[i] if len(s_data.get("T", [])) > i else "",
                    s_data.get("P", [])[i] if len(s_data.get("P", [])) > i else "",
                    s_data.get("flowrate", [])[i] if len(s_data.get("flowrate", [])) > i else "",
                ]
                for c in comps:
                    # Access the i-th value of the component's timeseries
                    comp_val = s_data.get("z", {}).get(c, [])[i] if len(s_data.get("z", {}).get(c, [])) > i else 0.0
                    row.append(comp_val)
                writer.writerow(row)


def save_results_json(results, fname="results.json"):
    """Save steady-state results to a JSON file."""
    os.makedirs("outputs", exist_ok=True)
    with open(os.path.join("outputs", fname), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=4)


def save_timeseries_json(results, fname="results_timeseries.json"):
    """
    Save dynamic results (time series) to a JSON file.
    results: dict {time: {stream: {...}}}
    """
    os.makedirs("outputs", exist_ok=True)
    with open(os.path.join("outputs", fname), "w", encoding="utf-8") as f:
        json.dump(results, f, indent=4)



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
    plt.savefig(os.path.join("outputs", "temps_" + fname))
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
    plt.savefig(os.path.join("outputs", "comps_" + fname))
    plt.close()



def plot_timeseries(results, fname="timeseries.png"):
    """Generate time-series plots for dynamic simulations."""
    streams = sorted(results.keys())
    if not streams:
        return
    
    # Assume all streams have the same time array from the first stream
    times = results[streams[0]].get("time", [])
    if not times:
        return

    # Temperature vs time
    plt.figure(figsize=(8, 5))
    for s_name in streams:
        if "T" in results[s_name]:
            plt.plot(times, results[s_name]["T"], label=s_name)
    plt.xlabel("Time [s]")
    plt.ylabel("Temperature [K]")
    plt.title("Stream Temperatures vs Time")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join("outputs", "temps_" + fname))
    plt.close()

    # Composition vs time (one plot per stream)
    comps = set()
    for s_data in results.values():
        if 'z' in s_data and isinstance(s_data['z'], dict):
            comps.update(s_data['z'].keys())
    comps = sorted(list(comps))

    for s_name in streams:
        if "z" not in results[s_name]:
            continue
        plt.figure(figsize=(8, 5))
        for comp in comps:
            if comp in results[s_name]["z"]:
                plt.plot(times, results[s_name]["z"][comp], label=comp)
        plt.xlabel("Time [s]")
        plt.ylabel("Mole Fraction")
        plt.title(f"Compositions vs Time ({s_name})")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join("outputs", f"comps_{s_name}_" + fname))
        plt.close()