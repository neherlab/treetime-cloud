#!/usr/bin/env bash

# Entrypoint for the base worker container dev mode.
# This will only prepare  virtual environement for the main worker containers.


set -o errexit
set -o nounset
set -o pipefail

cd /code/worker

mkdir -p ".cache"

# Tell main containers that the virtual environement is not ready
rm -f ".cache/venv.ready"

# Prepare virtual environment and install dependencies
yarn packages:install

# Tell main containers that the virtual environement is ready
touch ".cache/venv.ready"

# Important: this service is ephemerial, it only runs when the environement
# needs to be updated
exit 0
