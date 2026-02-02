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
    Test the particle sampling initialization

    Check if the number of fields and entries match and have valid results.

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared
#include <catch2/catch_template_test_macros.hpp>

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "basicItoPopeCloud.H"
#include "basicAerosolReactingPopeCloud.H"
#include "mmcStatusMessage.H"
#include "mmcVarSet.H"
#include "OFstream.H"
#include "IFstream.H"
#include "unitMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace particleStatisticsTest
{
    // Number of stored field entries for each particle type
    List<label> nFieldEntries = {1};
}

TEMPLATE_TEST_CASE("particleStatistics particle Test","[mmcClouds]",basicAerosolReactingPopeCloud)
{
    using particleType  = typename TestType::particleType;
    
    // Create the time and mesh
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object
    fvMesh mesh
    (
        IOobject
        (
            Foam::fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // Generate a cloud with particles to sample
    #include "createFieldsForMMC.H"

    INFO("Create the MMC cloud");
    TestType mmcCloud
    (
        "mmcCloud",
        mesh,
        U,
        DEff,
        rho,
        gradRhoDEff,
        Xi,
        false,
        false
    );


    // Add a single particle 
    label icell = -1;
    label itetface = -1;
    label itetpt = -1;
        
    mesh.findCellFacePt(mesh.C()[0],icell,itetface,itetpt);

    // Create a particle
    particleType p
    (
        mesh,
        Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
        icell,
        itetface,
        itetpt 
    );

    mmcCloud.setParticleProperties
    (
        p,      // particle
        1.0,    // mass
        1.0,    // weight
        0,      // patchi
        0,      // facei
        false   // init at boundary
    );

    wordList fields;
    auto varNames = p.getStatisticalDataNames(fields);
    auto varValues = p.getStatisticalData(fields);

    REQUIRE(varNames.size() == varValues.size());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

