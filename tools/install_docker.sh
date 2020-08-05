#!/usr/bin/env bash

export DEBIAN_FRONTEND=noninteractive

GOPATH="${GOPATH:-${HOME}/go}"

DOCKER_COMPOSE_VERSION=1.25.0
BUILDKIT_VERSION=0.7.1

CTOP_VERSION=0.7.3
DIVE_VERSION=0.9.2
LAZYDOCKER_VERSION=0.9
DRY_VERSION=0.10-beta.1

# Uninstall old versions:
sudo apt-get update -qq --yes
sudo apt-get remove -qq --yes \
containerd \
docker \
docker-engine \
docker.io \
runc \

sudo apt-get install -qq --install-recommends --yes \
apt-transport-https \
ca-certificates \
curl \
gnupg-agent \
software-properties-common \

# Add Docker’s official GPG key:
curl -fsSL "https://download.docker.com/linux/ubuntu/gpg" | sudo apt-key add -

# Add repo:
sudo sh -c "echo \"deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable\" > \"/etc/apt/sources.list.d/docker.list\""

# Install docker community edition:
sudo apt-get update -qq --yes
sudo apt-get install -qq --install-recommends --yes \
  build-essential \
  debhelper \
  devscripts \
  dh-make \
  docker-ce \
  docker-ce-cli \
  containerd.io \
  git \

# Allow current user to run docker without sudo:
sudo groupadd docker 2> /dev/null
sudo gpasswd -a ${USER} docker > /dev/null

sudo service docker restart

# docker-compose
pushd "/tmp" > /dev/null
  curl -fsSL "https://github.com/docker/compose/releases/download/${DOCKER_COMPOSE_VERSION}/docker-compose-$(uname -s)-$(uname -m)" -o "docker-compose"
  sudo chmod +x "docker-compose"
  sudo mv docker-compose /usr/local/bin/
popd > /dev/null

# buildkit
pushd "/tmp" > /dev/null
  curl -fsSL "https://github.com/moby/buildkit/releases/download/v${BUILDKIT_VERSION}/buildkit-v${BUILDKIT_VERSION}.linux-amd64.tar.gz" | tar -xz "bin/buildctl"
  sudo chmod +x "bin/buildctl"
  sudo mv bin/buildctl /usr/local/bin/
popd > /dev/null

# ctop
go get -u "github.com/ivan-aksamentov/ctop" >/dev/null
pushd "$GOPATH/src/github.com/ivan-aksamentov/ctop" >/dev/null
  make build >/dev/null
  sudo cp ctop /usr/local/bin/
popd > /dev/null

# dive
pushd "/tmp" > /dev/null
  curl -fsSL "https://github.com/wagoodman/dive/releases/download/v${DIVE_VERSION}/dive_${DIVE_VERSION}_linux_amd64.tar.gz" | tar -xz "dive"
  sudo chmod +x "dive"
  sudo mv dive /usr/local/bin/
popd > /dev/null

# lazydocker
pushd "/tmp" > /dev/null
  curl -fsSL "https://github.com/jesseduffield/lazydocker/releases/download/v${LAZYDOCKER_VERSION}/lazydocker_${LAZYDOCKER_VERSION}_Linux_x86_64.tar.gz" | tar -xz "lazydocker"
  chmod +x "lazydocker"
  sudo mv lazydocker /usr/local/bin/
popd > /dev/null

# dry
pushd "/tmp" > /dev/null
  curl -fsSL "https://github.com/moncho/dry/releases/download/v${DRY_VERSION}/dry-linux-amd64" -o "dry"
  sudo chmod +x "dry"
  sudo mv dry /usr/local/bin/
popd > /dev/null

# Reload group settings:
newgrp docker
