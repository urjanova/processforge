import csv
import json
import os

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
from loguru import logger

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


def generate_validation_excel(data_source, output_filename):
    """
    Generate a multi-sheet Excel validation report from simulation results.

    data_source: path to a CSV file or a pandas DataFrame.
    output_filename: path for the output .xlsx file.
    """
    # 1. Load Data
    if isinstance(data_source, str):
        df = pd.read_csv(data_source)
    else:
        df = data_source.copy()

    # Normalise the stream column name so both timeseries ('stream')
    # and steady-state ('Stream') CSVs are handled.
    if 'Stream' in df.columns and 'stream' not in df.columns:
        df.rename(columns={'Stream': 'stream'}, inplace=True)

    # Detect component columns (everything that is not a known metadata col)
    known_cols = {'time', 'stream', 'T [K]', 'P [Pa]', 'Phase',
                  'VaporFrac', 'flowrate'}
    comp_cols = [c for c in df.columns if c not in known_cols]

    # 2. Composition Integrity Check (Mass Balance)
    numeric_comp = df[comp_cols].apply(pd.to_numeric, errors='coerce')
    df['Total_Fraction'] = numeric_comp.sum(axis=1)
    df['Composition_Alert'] = np.where(
        np.isclose(df['Total_Fraction'], 1.0, atol=1e-5),
        "OK", "MASS LEAK"
    )
    mass_ok = (df['Composition_Alert'] == "OK").all()

    # 3. Unit Operation Performance (Pump Analysis)
    pump_check = pd.DataFrame()
    pump_ok = True
    temp_ok = True

    if 'stream' in df.columns and 'time' in df.columns:
        stream_names = df['stream'].unique()
        # Find pump pairs: streams named *before_pump* / *after_pump*
        pump_ins = sorted([s for s in stream_names if 'before_pump' in str(s)])
        pump_outs = sorted([s for s in stream_names if 'after_pump' in str(s)])

        for p_in, p_out in zip(pump_ins, pump_outs):
            df_in = df[df['stream'] == p_in].set_index('time')
            df_out = df[df['stream'] == p_out].set_index('time')
            common_idx = df_in.index.intersection(df_out.index)
            if common_idx.empty:
                continue
            pc = pd.DataFrame(index=common_idx)
            pc['Pump'] = f"{p_in} -> {p_out}"
            pc['Pressure_Gain_Pa'] = (
                df_out.loc[common_idx, 'P [Pa]'].values
                - df_in.loc[common_idx, 'P [Pa]'].values
            )
            pc['Temp_Rise_K'] = (
                df_out.loc[common_idx, 'T [K]'].values
                - df_in.loc[common_idx, 'T [K]'].values
            )
            pc['Pump_Status'] = np.where(
                pc['Pressure_Gain_Pa'] > 0, "Functional", "Broken"
            )
            pump_check = pd.concat([pump_check, pc])

        if not pump_check.empty:
            pump_ok = (pump_check['Pump_Status'] == "Functional").all()
            temp_ok = (pump_check['Temp_Rise_K'] >= 0).all()

    # 4. Build Executive Summary
    summary_rows = [
        {
            "Physical Law": "Conservation of Mass",
            "Logic": "Do chemical fractions add to 1.0?",
            "Status": "PASS" if mass_ok else "FAIL",
        }
    ]
    if not pump_check.empty:
        summary_rows.append({
            "Physical Law": "Pump Work (Pressure)",
            "Logic": "Does the pump increase pressure?",
            "Status": "PASS" if pump_ok else "FAIL",
        })
        summary_rows.append({
            "Physical Law": "Thermal Direction",
            "Logic": "Is the outlet temperature >= inlet?",
            "Status": "PASS" if temp_ok else "WARNING",
        })
    summary_df = pd.DataFrame(summary_rows)

    # 5. Export to Multi-Sheet Excel
    os.makedirs(os.path.dirname(output_filename) or ".", exist_ok=True)
    with pd.ExcelWriter(output_filename, engine='openpyxl') as writer:
        summary_df.to_excel(writer, sheet_name='1_EXECUTIVE_SUMMARY', index=False)
        if not pump_check.empty:
            pump_check.to_excel(writer, sheet_name='2_PUMP_PERFORMANCE')
        df.to_excel(writer, sheet_name='3_RAW_DATA_CHECKED', index=False)

    logger.info(f"Validation Report Generated: {output_filename}")