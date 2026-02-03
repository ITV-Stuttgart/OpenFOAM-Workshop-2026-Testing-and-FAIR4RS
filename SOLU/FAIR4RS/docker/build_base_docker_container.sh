#!/bin/bash

docker build                            \
       --rm                             \
       --build-arg UID=$(id -u)         \
       --build-arg GID=$(id -g)         \
       -t openfoam:v2406                \
       -f docker/Dockerfile.createOFImage \
       . 

