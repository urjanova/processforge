#!/bin/bash

set -euo pipefail

usage() {
	cat <<'EOF'
Usage:
  ./build_and_publish.sh test      # Build + upload to TestPyPI only
  ./build_and_publish.sh release   # Build + upload to PyPI only
  ./build_and_publish.sh both      # Build + upload to TestPyPI, then PyPI

Scripts:
  ./build_and_upload_testpypi.sh
  ./build_and_upload_pypi.sh
EOF
}

mode="${1:-}"

case "$mode" in
	test)
		./build_and_upload_testpypi.sh
		;;
	release)
		./build_and_upload_pypi.sh
		;;
	both)
		./build_and_upload_testpypi.sh
		./build_and_upload_pypi.sh
		;;
	*)
		usage
		exit 1
		;;
esac