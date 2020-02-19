#!/usr/bin/env bash

# Create directories necessary for docker operation

set -o errexit
set -o nounset
set -o pipefail

mkdir -p \
  .volumes/filestore \
  .volumes/nginx/logs \
  .volumes/nginx/tmp/client_body_temp_path \
  .volumes/nginx/tmp/proxy_temp_path \
  .volumes/nginx/cache \
  .volumes/taskqueue/config/generated \
  .volumes/taskqueue/home \
  .volumes/taskqueue/log \
  .volumes/taskqueue/schema
