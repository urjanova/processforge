#!/usr/bin/env bash
# Railway mounts its Volume at /data, owned by root. The application runs as
# the unprivileged MAMBA_USER; run artifacts are written to the scratch dir
# set by PROCESSFORGE_OUTPUT_DIR (the Docker image defaults to /tmp/processforge,
# so scratch stays container-local and ephemeral). Ensure /data is writable by
# that user (used for cross-section data), then drop privileges and exec the
# server/command.
set -euo pipefail

DATA_ROOT="${OPENMC_DATA_ROOT:-/data}"

if [ "$(id -u)" = "0" ] && [ -n "${MAMBA_USER:-}" ]; then
    # Make the mounted volume writable by the runtime user. Harmless when no
    # volume is mounted (the anonymous VOLUME dir is root-owned and empty).
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
