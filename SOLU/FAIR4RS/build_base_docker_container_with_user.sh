#!/bin/bash

docker build                            \
       --rm                             \
       --build-arg UID=$(id -u)         \
       --build-arg GID=$(id -g)         \
       -t openfoam:v2406user            \
       -f docker/Dockerfile.createOFImage_with_user \
       . 

