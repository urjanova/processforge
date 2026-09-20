#!/bin/bash
set -euo pipefail

# Geant4 does not require a downloadable cross-section data archive like OpenMC.
# This entrypoint simply forwards to the provider API server.
exec "$@"