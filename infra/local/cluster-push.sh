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

REGISTRY_ADDR=$(minikube --profile="${CLUSTER_NAME}" ip)
REGISTRY_PORT=5000
REGISTRY="${REGISTRY_ADDR}:${REGISTRY_PORT}"

images=(
  "treetime-prod-api"
  "treetime-prod-reverseproxy"
  "treetime-prod-taskqueue"
  "treetime-prod-web"
  "treetime-prod-worker"
)


# print_troubleshooting() {
# cat << EOF

# ---

# If you keep receiving an error "http: server gave HTTP response to HTTPS client",
# try to add the following entry to your '~/.docker/daemon.json' file:

# {
#     "insecure-registries" : [ "${REGISTRY}" ]
# }

# and restart your docker service, e.g.:

# $ sudo systemctl restart docker

# See: https://github.com/docker/distribution/issues/1874#issuecomment-237194314

# EOF
# }

for image in ${images[@]}; do
  NEW_TAG="${REGISTRY}/${image}:latest"
  # docker tag "${image}" "${NEW_TAG}"
  minikube cache add "${image}"

  # if ! docker push "${NEW_TAG}" ;then
  #   print_troubleshooting
  #   exit 1
  # fi
done

minikube cache reload
