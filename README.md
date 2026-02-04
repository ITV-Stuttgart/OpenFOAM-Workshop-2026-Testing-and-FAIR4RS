# OpenFOAM Workshop 2026

Ensuring the reliability, reproducibility, and long-term sustainability of
scientific software is becoming increasingly essential as computational
research grows in complexity and scale. OpenFOAM, as one of the most widely
used open-source CFD frameworks, provides immense flexibility, but this
flexibility also increases the need for robust testing practices and
transparent, reproducible research workflows. This workshop, “Testing
Strategies for OpenFOAM and FAIR4RS Principles,” introduces participants to
modern methodologies for improving software quality and reproducibility in the
OpenFOAM ecosystem.

The first part of the workshop focuses on testing strategies, following
concepts outlined in the paper [Testing Strategies for OpenFOAM](https://doi.org/10.51560/ofj.v5.134) 
of Gärtner et al.. Participants will learn why testing and testable code are
critical for maintaining correctness, facilitating extension, and supporting
collaborative development. We will introduce different categories of tests --
static tests, unit tests, regression tests, and integration tests -- and define
each within the specific context of OpenFOAM’s architecture and solver
development. After a conceptual introduction, participants will gain hands-on
experience integrating testing frameworks into an example solver, including
configuring static tests in git and setting up a unit-testing workflow using
Catch2. By the end of this module, attendees will understand how systematic
testing increases confidence in new models, reduces maintenance burden, and
enables more robust community contributions.

The second part of the workshop introduces the FAIR Principles for Research
Software (FAIR4RS), emphasizing their relevance for computational research and
community-driven software such as OpenFOAM. We will discuss what it means for
research software to be Findable, Accessible, Interoperable, and Reusable, and
highlight practical mechanisms to support these goals. A major emphasis will be
placed on reproducibility: participants will explore different strategies for
capturing and recreating computational environments, including Docker and
related containerisation tools, while also discussing their limitations—for
example, challenges associated with multi-node or HPC deployments. The workshop
will then demonstrate practical steps to enhance reproducibility in OpenFOAM
workflows, such as embedding git version information into solver log files,
managing dependencies, and documenting computational setups to enable long-term
replicability.


## Contents

This directory contains the case files for the OpenFOAM Workshop 2026 at the 
TU Eindhoven. 
The workshop imparts knowledge about testing OpenFOAM applications with the 
unit test framework [Catch2](https://github.com/catchorg/Catch2), presenting a 
test case for the perfect gas equation of state and a coupled test for an
evaporating Lagrangian particle with a custom evaporation function. Further, a 
case study with the containerization tool docker and apptainer is given, to
showcast how a custom solver can be packaged for sharing, archiving, and 
reproducing results. 

## License
Attribution-NonCommercial-ShareAlike 4.0 International


