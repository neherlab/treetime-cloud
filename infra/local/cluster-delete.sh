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

if minikube status --profile="${CLUSTER_NAME}" >/dev/null; then
  minikube delete --profile "${CLUSTER_NAME}"
else
  echo "Cluster '${CLUSTER_NAME}' does not exist"
  exit 1
fi
