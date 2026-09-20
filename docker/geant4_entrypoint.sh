#!/bin/bash
set -euo pipefail

# Allow Geant4 to be served locally or remotely
# The provider listens on port 9003

# Initialize the Geant4 provider and start serving
exec micromamba run -n base python -m processforge.api.serve "$@"