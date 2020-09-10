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

kind delete cluster --name "${CLUSTER_NAME}"
