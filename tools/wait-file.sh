#!/usr/bin/env bash

# Blocks execution, waiting for a given a file to be created or updated.

set -o errexit
set -o nounset
set -o pipefail

# Path to watch for
FILEPATH="${1:-}"

if [ -z "${FILEPATH}" ]; then
  echo "Usage: ${0} <filepath>"
  exit 1
fi

FILENAME=$(basename ${FILEPATH})
DIRNAME=$(realpath $(dirname "${PWD}/${FILEPATH}"))

if [ ! -f "${FILEPATH}" ]
then
  echo "Waiting for ${DIRNAME}/${FILENAME}"
  while read filename; do
    if [ "$filename" = ${FILENAME} ]; then
      break
    fi
  done < <(inotifywait -e create,open --format '%f' --quiet ${DIRNAME} --monitor)
fi
