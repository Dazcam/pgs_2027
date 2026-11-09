#!/usr/bin/env bash
#--------------------------------------------------------------------------------------
#
#    Pull Singularity / Apptainer containers
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  00SETUP | Rule: get_containers
#
# Purpose:   Pull each container to its target .sif path; skip any that already exist
#
# Usage:     setup_pull_containers.sh <target.sif> <uri> [<target.sif> <uri> ...]
#            URIs can be oras:// or docker:// (docker images are converted to .sif)
#
#--------------------------------------------------------------------------------------

set -euo pipefail

if [ $# -eq 0 ] || [ $(( $# % 2 )) -ne 0 ]; then
    echo "ERROR: arguments must be <target.sif> <uri> pairs"
    exit 1
fi

if command -v apptainer &> /dev/null; then
    CMD="apptainer"
elif command -v singularity &> /dev/null; then
    CMD="singularity"
else
    echo "ERROR: neither apptainer nor singularity found"
    exit 1
fi

echo "Using: $CMD"

while [ $# -gt 0 ]; do
    target="$1"
    uri="$2"
    shift 2

    if [ -f "$target" ]; then
        echo "SKIP: $target already exists"
        continue
    fi

    mkdir -p "$(dirname "$target")"
    echo "Pulling: $uri -> $target"
    "$CMD" pull "$target" "$uri"
    ls -lh "$target"
done

echo "Done."

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
