#!/bin/bash

# Source OpenFOAM
. $foamDotFile

CASE_DIR="${CASE_DIR:-/home/${OF_USER}/workspace}"

mkdir -p $CASE_DIR

cd $CASE_DIR

if [[ "$#" -gt 0 ]]; then
    exec dropletEvapFoam "$@"
else
    exec dropletEvapFoam
fi

