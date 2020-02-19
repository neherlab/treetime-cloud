#!/usr/bin/env bash

# Create directories necessary for docker operation

set -o errexit
set -o nounset
set -o pipefail

mkdir -p \
  .volumes/filestore \
  .volumes/taskqueue/config/generated \
  .volumes/taskqueue/home \
  .volumes/taskqueue/log \
  .volumes/taskqueue/schema
