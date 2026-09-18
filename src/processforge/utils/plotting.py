"""Terminal plotting helpers for processforge run results."""

from __future__ import annotations

import numpy as np
import plotext as plt

from ..result import _convert_value, _scalar_from_sequence

# Default terminal plot dimensions.
_DEFAULT_WIDTH = 80
_DEFAULT_HEIGHT = 18
_MAX_WIDTH = 120
_MAX_HEIGHT = 40


def _plotext_figure(width: int = _DEFAULT_WIDTH, height: int = _DEFAULT_HEIGHT):
    """Return a fresh plotext figure with a sensible terminal size."""
    fig = plt.figure
    fig.clear()
    fig.plot_size(width, height)
    fig.theme("dark")
    return fig


def _shorten_labels(labels, max_len: int = 14):
    """Shorten long labels while preserving readability."""
    return [lab if len(lab) <= max_len else lab[: max_len - 1] + "…" for lab in labels]


def _bar_figure(n_categories: int, max_label_len: int):
    """Create a figure sized appropriately for a bar chart.

    Uses horizontal bars when there are many categories or long labels,
    vertical bars otherwise.
    """
    use_horizontal = n_categories > 6 or max_label_len > 12
    if use_horizontal:
        height = min(_MAX_HEIGHT, max(_DEFAULT_HEIGHT, n_categories * 2 + 6))
        return _plotext_figure(width=_DEFAULT_WIDTH, height=height), "horizontal"
    width = min(_MAX_WIDTH, max(_DEFAULT_WIDTH, n_categories * 12))
    return _plotext_figure(width=width, height=_DEFAULT_HEIGHT), "vertical"


def plot_results_terminal(results, title=""):
    """Render steady-state stream results as terminal bar charts."""
    streams = [k for k, v in results.items() if isinstance(v, dict) and "T" in v]
    if not streams:
        return

    temp_values = [_scalar_from_sequence(results[s].get("T")) or 0.0 for s in streams]
    labels = _shorten_labels(streams)
    max_len = max((len(lab) for lab in labels), default=0)

    fig, orientation = _bar_figure(len(streams), max_len)
    fig.draw(fig.bar(labels, temp_values, orientation=orientation, width=0.5))
    fig.title(title or "Stream Temperatures (Steady State)")
    if orientation == "horizontal":
        fig.label("Temperature (K)", axis="x")
        fig.label("Stream", axis="y")
    else:
        fig.label("Stream", axis="x")
        fig.label("Temperature (K)", axis="y")
    fig.show()
    fig.clear()

    comps = set()
    for s in streams:
        z = results[s].get("z")
        if isinstance(z, dict):
            comps.update(z.keys())
    comps = sorted(comps)
    if not comps:
        return

    values = []
    for comp in comps:
        values.append([
            _scalar_from_sequence(results[s].get("z", {}).get(comp, 0.0)) or 0.0
            for s in streams
        ])

    fig, orientation = _bar_figure(len(streams), max_len)
    fig.draw(fig.bar(labels, values, orientation=orientation, stacked=True, width=0.5))
    fig.title(title or "Stream Compositions (Steady State)")
    if orientation == "horizontal":
        fig.label("Mole fraction", axis="x")
        fig.label("Stream", axis="y")
    else:
        fig.label("Stream", axis="x")
        fig.label("Mole fraction", axis="y")
    fig.show()
    fig.clear()


def plot_timeseries_terminal(results, title=""):
    """Render dynamic stream results as terminal line plots."""
    streams = sorted(k for k, v in results.items() if isinstance(v, dict))
    if not streams:
        return

    times = results[streams[0]].get("time", [])
    if not times:
        return

    n_points = len(times)
    width = min(_MAX_WIDTH, max(_DEFAULT_WIDTH, n_points // 2))

    fig = _plotext_figure(width=width, height=_DEFAULT_HEIGHT)
    for s_name in streams:
        if "T" in results[s_name]:
            fig.draw(fig.signal(times, results[s_name]["T"]).label(s_name))
    fig.title(title or "Stream Temperatures vs Time")
    fig.label("Time (s)", axis="x")
    fig.label("Temperature (K)", axis="y")
    fig.show()
    fig.clear()

    comps = set()
    for s_data in results.values():
        if "z" in s_data and isinstance(s_data["z"], dict):
            comps.update(s_data["z"].keys())
    comps = sorted(comps)

    for s_name in streams:
        if "z" not in results[s_name]:
            continue
        fig = _plotext_figure(width=width, height=_DEFAULT_HEIGHT)
        for comp in comps:
            if comp in results[s_name]["z"]:
                fig.draw(fig.signal(times, results[s_name]["z"][comp]).label(comp))
        fig.title(f"Compositions vs Time ({s_name})")
        fig.label("Time (s)", axis="x")
        fig.label("Mole fraction", axis="y")
        fig.show()
        fig.clear()


def plot_zarr_summary_terminal(summary: dict, title: str = "") -> None:
    """Render a terminal visualization from a Zarr summary dict."""
    if not summary.get("present"):
        return

    mode = summary.get("mode", "steady")
    streams = summary.get("streams", {})

    if mode == "dynamic":
        # Reconstruct a results-like dict from Zarr summary for plotting.
        results: dict = {}
        for s_name, sdata in streams.items():
            fields = sdata.get("fields", {})
            results[s_name] = {k: v.get("value") for k, v in fields.items()}
        plot_timeseries_terminal(results, title=title)
    else:
        # Steady: show final timestep values.
        results = {}
        for s_name, sdata in streams.items():
            fields = sdata.get("fields", {})
            results[s_name] = {
                k: _scalar_from_sequence(v.get("value"))
                for k, v in fields.items()
            }
        plot_results_terminal(results, title=title)


def plot_results_to_terminal(results, title: str = "", mode: str = "steady") -> None:
    """Dispatch to the appropriate terminal plotter for run results."""
    if mode == "dynamic":
        plot_timeseries_terminal(results, title=title)
    else:
        plot_results_terminal(results, title=title)
