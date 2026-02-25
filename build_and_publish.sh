#!/bin/bash

# install build tools if you haven't already
pip install --upgrade build twine

# from the repo root
python -m build # creates dist/processforge‑X.Y.Z.tar.gz  dist/processforge‑X.Y.Z‑py3-none-any.whl

# optionally inspect the contents
ls dist

# if you have not registered on pypi.org, create an account there and
# configure ~/.pypirc or export TWINE_USERNAME/TWINE_PASSWORD

# upload to the test PyPI first (strongly recommended)
python -m twine upload --repository testpypi dist/*

# when you're happy, upload to the real PyPI
python -m twine upload dist/*

# you can also specify --repository pypi explicitly
# to upload only a single file, e.g. twine upload dist/processforge-0.2.0-py3-none-any.whl