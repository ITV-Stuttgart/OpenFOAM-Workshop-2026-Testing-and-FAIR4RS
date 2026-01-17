/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

Description
    Macros for creating basic particle thermo mixtures

\*---------------------------------------------------------------------------*/

#ifndef makeThermo_H
#define makeThermo_H

#include "particleThermoMixture.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie)          \
                                                                               \
    typedef                                                                    \
        Transport                                                              \
        <                                                                      \
            species::thermo                                                    \
            <                                                                  \
                Thermo                                                         \
                <                                                              \
                    EqnOfState                                                 \
                    <                                                          \
                        Specie                                                 \
                    >                                                          \
                >,                                                             \
                Type                                                           \
            >                                                                  \
        >                                                                      \
        Transport##Type##Thermo##EqnOfState##Specie;

#define defineParticleThermoMixture(ThermoPhys)                                \
                                                                               \
    typedef particleThermoMixture<ThermoPhys>                                  \
        ParticleThermoMixture##ThermoPhys;                                     \
                                                                               \
    defineTemplateTypeNameAndDebugWithName                                     \
    (                                                                          \
        ParticleThermoMixture##ThermoPhys,                                     \
        ("particleThermoMixture<" + ThermoPhys::typeName() + ">").c_str(),     \
        0                                                                      \
    );


#define addParticleThermoMixture(BaseParticleThermoMixture,ThermoPhys) \
                                                                               \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        BaseParticleThermoMixture,                                             \
        ParticleThermoMixture##ThermoPhys,                                     \
        dictionary                                                             \
    );                                                                         \


#define makeParticleThermoMixtures(Transport,Type,Thermo,EqnOfState,Specie)    \
    typedefThermoPhysics(Transport,Type,Thermo,EqnOfState,Specie)              \
                                                                               \
    defineParticleThermoMixture(Transport##Type##Thermo##EqnOfState##Specie)   \
                                                                               \
    addParticleThermoMixture(baseParticleThermoMixture,Transport##Type##Thermo##EqnOfState##Specie)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
