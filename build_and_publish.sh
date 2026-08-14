#!/bin/bash

set -euo pipefail

VERSION="$(grep -m1 '^version' pyproject.toml | sed -E 's/version *= *"([^"]+)".*/\1/')"
TAG="v${VERSION}"

tag_release() {
	if git rev-parse -q --verify "refs/tags/${TAG}" >/dev/null; then
		echo "Tag ${TAG} already exists; skipping."
	else
		echo "Creating tag ${TAG}..."
		git tag -a "${TAG}" -m "Release ${VERSION}"
		git push origin "${TAG}"
	fi
}

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
		tag_release
		;;
	release)
		./build_and_upload_pypi.sh
		tag_release
		;;
	both)
		./build_and_upload_testpypi.sh
		./build_and_upload_pypi.sh
		tag_release
		;;
	*)
		usage
		exit 1
		;;
esac