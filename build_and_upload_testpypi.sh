#!/bin/bash

set -euo pipefail

[ -d dist ] && rm -r dist/

uv build
ls -lh dist

: "${TEST_PYPI_API_TOKEN:?Set TEST_PYPI_API_TOKEN before running this script}"

uvx twine upload \
  --verbose \
  --non-interactive \
  --repository-url https://test.pypi.org/legacy/ \
  -u __token__ \
  -p "$TEST_PYPI_API_TOKEN" \
  dist/*
