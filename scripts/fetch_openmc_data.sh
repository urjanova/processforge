#!/usr/bin/env bash
# Download OpenMC data files to OPENMC_DATA_ROOT if not already present.
#
# Environment variables (all optional — script is a no-op if none are set):
#
#   OPENMC_DATA_ROOT   Local directory for cross sections + geometry.
#                      Defaults to /data.
#
#   OPENMC_DATA_URL    URL to a .tar.gz archive of the cross-section library
#                      (e.g. an ENDF/B-VIII.0 HDF5 tarball from the OpenMC
#                      data releases page). Downloaded once; skipped on
#                      subsequent container starts if cross_sections.xml
#                      already exists in OPENMC_DATA_ROOT.
#                      After download the script exports OPENMC_CROSS_SECTIONS
#                      automatically — do not set that variable manually.
#
#   OPENMC_DAGMC_URL   URL to a single DAGMC .h5m geometry file. Downloaded
#                      once to $OPENMC_DATA_ROOT/geometry/. The filename is
#                      taken from the URL basename.
#
# Docker usage:
#   docker run -v /host/openmc_data:/data processforge run flowsheet.json
#
# Railway usage:
#   Attach a Railway Volume at /data, then set OPENMC_DATA_URL and
#   OPENMC_DAGMC_URL in the Environment Variables tab. The script
#   downloads on first deploy; subsequent deploys skip the download.
set -euo pipefail

DATA_ROOT="${OPENMC_DATA_ROOT:-/data}"
mkdir -p "$DATA_ROOT"

# ---------------------------------------------------------------------------
# Cross sections
# ---------------------------------------------------------------------------
XS_FILE="$DATA_ROOT/cross_sections.xml"

if [ -n "${OPENMC_DATA_URL:-}" ] && [ ! -f "$XS_FILE" ]; then
    echo "[processforge] Downloading cross sections from $OPENMC_DATA_URL ..."
    curl -fsSL "$OPENMC_DATA_URL" -o /tmp/openmc_data.tar.gz
    tar -xzf /tmp/openmc_data.tar.gz -C "$DATA_ROOT" --strip-components=1
    rm -f /tmp/openmc_data.tar.gz
    echo "[processforge] Cross sections ready at $XS_FILE"
fi

# Export OPENMC_CROSS_SECTIONS so OpenMC can find the library at runtime.
# This works whether the data was just downloaded or was already present on the volume.
if [ -f "$XS_FILE" ]; then
    export OPENMC_CROSS_SECTIONS="$XS_FILE"
fi

# ---------------------------------------------------------------------------
# DAGMC geometry
# ---------------------------------------------------------------------------
if [ -n "${OPENMC_DAGMC_URL:-}" ]; then
    DAGMC_DIR="$DATA_ROOT/geometry"
    mkdir -p "$DAGMC_DIR"
    # Strip query string to get a clean filename
    DAGMC_FILE="$DAGMC_DIR/$(basename "${OPENMC_DAGMC_URL%%\?*}")"
    if [ ! -f "$DAGMC_FILE" ]; then
        echo "[processforge] Downloading DAGMC geometry from $OPENMC_DAGMC_URL ..."
        curl -fsSL "$OPENMC_DAGMC_URL" -o "$DAGMC_FILE"
        echo "[processforge] Geometry ready at $DAGMC_FILE"
    fi
fi
