#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

set -x

# Gets directory where this script is located
THIS_DIR=$(cd $(dirname "${BASH_SOURCE[0]}"); pwd)

# Loads variables from file
if [ -f "${THIS_DIR}/settings" ]; then
  export $(echo $(cat "${THIS_DIR}/settings" | sed 's/#.*//g'| xargs))
fi

# Creates a cluster network, if not already present
maybe-create-cluster-network() {
  if ! docker network ls --format "{{.Name}}" | grep -w -q "${CLUSTER_NETWORK_NAME}"; then
    docker network create "${CLUSTER_NETWORK_NAME}" --subnet="${CLUSTER_NETWORK_SUBNET}"
  fi
}

# Runs a registry container, if not already running
maybe-create-registry() {
  RUNNING=$(docker inspect -f '{{.State.Running}}' docker-registry 2>/dev/null || true )
  if [ "${RUNNING}" != "true" ]; then
    docker run -d --restart=always -p "${REGISTRY_PORT_EXTERNAL}:${REGISTRY_PORT_INTERNAL}" --name "${REGISTRY_NAME}" registry:2
  fi
}

# Lists containers connected to the cluster network and connects registry container if it's not on the list
maybe-connect-registry-to-cluster-network() {
  CONTAINERS=$(docker network inspect "${CLUSTER_NETWORK_NAME}" -f "{{range .Containers}}{{.Name}} {{end}}" | xargs echo -n)
  if [[ ! " ${CONTAINERS[@]} " =~ " ${REGISTRY_NAME} " ]]; then
    docker network connect "${CLUSTER_NETWORK_NAME}" "${REGISTRY_NAME}"
  fi
}

# Creates a cluster using 'kind', if not already present.
# Uses config file with env variable substitutions.
maybe-create-cluster() {
  if ! kind get clusters -q | grep -w -q "${CLUSTER_NAME}"; then
    envsubst < "${THIS_DIR}/kind.yml" | kind create cluster --name "${CLUSTER_NAME}" --config=-
  fi
}

maybe-create-cluster-network
maybe-create-registry
maybe-connect-registry-to-cluster-network
maybe-create-cluster
