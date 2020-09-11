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

if ! command -v docker version &> /dev/null; then
  >&2 printf "Error: docker not found.\nFollow the installation guide at:\nhttps://docs.docker.com/get-docker"
  exit 1
fi

if ! command -v kubectl version &> /dev/null; then
  >&2 printf "Error: kubectl not found.\nFollow the installation guide at:\nhttps://kubernetes.io/docs/tasks/tools/install-kubectl/"
  exit 1
fi

if ! command -v minikube version &> /dev/null; then
  >&2 printf "Error: minikube not found.\nFollow the installation guide at:\nhttps://minikube.sigs.k8s.io/docs/start/"
  exit 1
fi

export KUBECTL_VERSION=$(kubectl version --client --short | awk -F ' ' '{print $3}' | xargs echo -n)

minikube start \
--profile="${CLUSTER_NAME}" \
--kubernetes-version="${KUBECTL_VERSION}" \
--addons=dashboard \
--cpus=${CLUSTER_NUM_CPUS} \
--memory=${CLUSTER_MEMORY} \
--nodes=${CLUSTER_NUM_NODES} \
--wait=all \
--interactive=false

# --addons=default-storageclass \
# --addons=helm-tiller \
# --addons=ingress \
# --addons=logviewer \
# --addons=registry \
# --addons=storage-provisioner \
