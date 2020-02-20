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
cat <<EOT >> ".cache/venv.ready"
This file is created by '$(pwd)/${0}'

Existence of this feale means that the python-poetry virtual environment is
ready and services waiting for it can proceed.

Here is some additional information about the virtual environment:
$(poetry env info)
EOT

# Important: this service is ephemeral, it only runs when the environement
# needs to be updated, updates it and the exits.
exit 0
