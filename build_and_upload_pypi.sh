#!/bin/bash

set -euo pipefail

[ -d dist ] && rm -r dist/

python -m pip install --upgrade build twine
python -m build
ls -lh dist

: "${PYPI_API_TOKEN:?Set PYPI_API_TOKEN before running this script}"

python -m twine upload \
  --verbose \
  --non-interactive \
  --repository-url https://upload.pypi.org/legacy/ \
  -u __token__ \
  -p "$PYPI_API_TOKEN" \
  dist/*
