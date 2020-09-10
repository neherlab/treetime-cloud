#!/usr/bin/env bash

# Create directories necessary for docker operation

set -o errexit
set -o nounset
set -o pipefail

mkdir -p \
  .volumes/development/filestore \
  .volumes/development/nginx/logs \
  .volumes/development/nginx/tmp/client_body_temp_path \
  .volumes/development/nginx/tmp/proxy_temp_path \
  .volumes/development/nginx/cache \
  .volumes/development/taskqueue/config/generated \
  .volumes/development/taskqueue/home \
  .volumes/development/taskqueue/log \
  .volumes/development/taskqueue/schema
