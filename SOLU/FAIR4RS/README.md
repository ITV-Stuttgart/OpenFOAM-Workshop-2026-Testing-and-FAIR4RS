# Reproducibility and Archiving OpenFOAM with Docker and Apptainer

> [!NOTE] 
> This file is part of the OpenFOAM Workshop 2026

Based on the previous testing steps using the Catch2 example, this case
includes a Lagrangian library and an OpenFOAM solver that uses this library,
called `dropletEvapFoam`. In the following, this custom solver is packaged
using Docker and Apptainer to create container images that allow sharing and
archiving a specific code version for reproducible results.

---

# Docker

Docker is a widely used container technology and therefore a good choice for
packaging custom OpenFOAM solvers in a distributable form. The following
sections describe the steps required to create a Docker container containing
OpenFOAM together with a custom solver.

## Create a Docker image

When creating a Docker image with OpenFOAM and a custom library or solver, it
is recommended to split the process into two steps. First, create a base Docker
image with OpenFOAM installed for the operating system of your choice. In a
second step, use this base image to build the final image that includes
OpenFOAM and your custom solver.

Docker images are created using a so-called `Dockerfile`. This file contains
the instructions required to install OpenFOAM and compile your library or
solver, similar to a bash script. Further details on Dockerfiles can be found
here: https://docs.docker.com/reference/dockerfile/

### Create the Docker image for this workshop

For this workshop, the first and most time-consuming step has already been
completed. A ready-to-use Docker image with OpenFOAM v2406 installed is
available at: jgaertner2/openfoam:v2406

To install and package the custom solver and library, run the provided script
`build_dropletEvapFoam_docker_container.sh`. This script pulls the OpenFOAM
v2406 Docker image, copies the custom solver and library into the container
using the Docker `COPY` command, and compiles them using `wmake`.

Finally, an entrypoint script is included so that the solver is executed
automatically in the case directory when the container is run.

## Run the Docker image

To run the Docker image with the installed solver, the OpenFOAM case directory
must be mounted into the container:

```bash
docker run -v ${PWD}:/home/ofuser/workspace dropletevapfoam:v1.0
```

This command assumes that it is executed from within the case directory.


# Apptainer

Docker containers offer many advantages, but they are typically unsuitable for
HPC systems due to their requirement for superuser privileges. Apptainer is an
alternative container solution designed for HPC environments and can directly
use Docker images.

## Generate an Apptainer image

Apptainer creates a single SIF image that contains everything required to run
the application. This makes it well suited for archiving solvers for long-term
reproducibility, especially when dependencies or compilation environments may
no longer be available.

To generate the Apptainer image, run:
```bash
apptainer build apptainer/solver.sif apptainer/dropletEvapFoam.def
```

On Ubuntu 22.04 and earlier, Apptainer requires superuser privileges. Alternatively, it can be used with the --fakeroot option:
```bash
apptainer build --fakeroot apptainer/solver.sif apptainer/dropletEvapFoam.def
```

## Run the Apptainer Image

To execute the Apptainer image, navigate to the OpenFOAM case directory and run:
```bash
apptainer run --fakeroot apptainer/solver.sif
```
This will start the solver using the packaged environment.
Without fakeroot the files written will still have root ownership.




