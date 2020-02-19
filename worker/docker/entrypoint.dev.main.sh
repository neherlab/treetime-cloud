#!/usr/bin/env bash

# Entrypoint for the main worker container dev mode.
# This will run the actual worker application.
# Base container is expected to run before that.

set -o errexit
set -o nounset
set -o pipefail

cd /code/worker

# Wait until base container prepares the virtual environement
/tools/wait-file.sh ".cache/venv.ready"

# Wait until taskqueue service is available
dockerize -wait "tcp://treetime-dev-taskqueue:5672" -timeout 60s >& /dev/null

# Run the application
yarn dev
