#!/bin/bash

docker build                            \
       --rm                             \
       -t dropletevapfoam:v1.0          \
       -f docker/Dockerfile.dropletEvapFoam    \
       .

