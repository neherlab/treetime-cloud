#!/usr/bin/env bash

# Create directories necessary for docker operation

set -o errexit
set -o nounset
set -o pipefail

mkdir -p \
  .volumes/production/filestore \
  .volumes/production/nginx/logs \
  .volumes/production/nginx/tmp/client_body_temp_path \
  .volumes/production/nginx/tmp/proxy_temp_path \
  .volumes/production/nginx/cache \
  .volumes/production/taskqueue/config/generated \
  .volumes/production/taskqueue/home \
  .volumes/production/taskqueue/log \
  .volumes/production/taskqueue/schema

# Prevent rabbitmq from crashing with error ".erlang.cookie must be accessible by owner only"
chmod 600 services/taskqueue/prod/var/lib/rabbitmq/.erlang.cookie
