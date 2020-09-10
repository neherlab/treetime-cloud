#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

# Gets directory where this script is located
THIS_DIR=$(cd $(dirname "${BASH_SOURCE[0]}"); pwd)

# Loads variables from file
if [ -f "${THIS_DIR}/settings" ]; then
  export $(echo $(cat "${THIS_DIR}/settings" | sed 's/#.*//g'| xargs))
fi

REGISTRY=localhost:${REGISTRY_PORT_EXTERNAL}

images=(
  "treetime-prod-api"
  "treetime-prod-reverseproxy"
  "treetime-prod-taskqueue"
  "treetime-prod-web"
  "treetime-prod-worker"
)

for image in ${images[@]}; do
  NEW_TAG="${REGISTRY}/${image}:latest"
  docker tag "${image}" "${NEW_TAG}"
  docker push "${NEW_TAG}"
done
