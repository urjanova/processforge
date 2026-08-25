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
# Docker usage:
#   docker run -v /host/openmc_data:/data processforge run flowsheet.json
#
# Railway usage:
#   Attach a Railway Volume at /data, then set OPENMC_DATA_URL
#   in the Environment Variables tab. The script
#   downloads on first deploy; subsequent deploys skip the download.
set -euo pipefail

DATA_ROOT="${OPENMC_DATA_ROOT:-/data}"
mkdir -p "$DATA_ROOT"

# ---------------------------------------------------------------------------
# Cross sections
# ---------------------------------------------------------------------------
XS_DIR="$DATA_ROOT/cross_sections"
XS_FILE="$XS_DIR/cross_sections.xml"

if [ -n "${OPENMC_DATA_URL:-}" ] && [ ! -f "$XS_FILE" ]; then
    echo "[processforge] Downloading cross sections from $OPENMC_DATA_URL ..."
    curl -fsSL "$OPENMC_DATA_URL" -o /tmp/openmc_data.tar.gz
    mkdir -p "$XS_DIR"
    tar -xzf /tmp/openmc_data.tar.gz -C "$XS_DIR" --strip-components=1
    rm -f /tmp/openmc_data.tar.gz
    echo "[processforge] Cross sections ready at $XS_FILE"
fi

# Export OPENMC_CROSS_SECTIONS so OpenMC can find the library at runtime.
# This works whether the data was just downloaded or was already present on the volume.
if [ -f "$XS_FILE" ]; then
    export OPENMC_CROSS_SECTIONS="$XS_FILE"
fi

# ---------------------------------------------------------------------------
# Drop privileges for the long-running server/command
# ---------------------------------------------------------------------------
# The data-prep work above (mkdir /data, curl, tar) needs write access to the
# mounted /data volume, which is owned by root on Railway. We therefore run the
# entrypoint as root. The server itself only needs to *read* the cross sections
# and writes its artifacts to PROCESSFORGE_OUTPUT_DIR (/tmp/processforge), so we
# drop back to the unprivileged MAMBA_USER before exec'ing the command.
if [ "$(id -u)" = "0" ] && [ -n "${MAMBA_USER:-}" ]; then
    # Make the mounted volume writable by the runtime user. The server only
    # needs to *read* cross sections, but chowning the root lets it also write
    # outputs there if PROCESSFORGE_OUTPUT_DIR is pointed at /data.
    if [ -d "$DATA_ROOT" ]; then
        chown "$MAMBA_USER" "$DATA_ROOT" 2>/dev/null || true
    fi
    if command -v runuser >/dev/null 2>&1; then
        exec runuser -u "$MAMBA_USER" -- "$@"
    elif command -v su >/dev/null 2>&1; then
        exec su -s /bin/bash "$MAMBA_USER" -c 'exec "$@"' bash "$@"
    fi
fi

exec "$@"
