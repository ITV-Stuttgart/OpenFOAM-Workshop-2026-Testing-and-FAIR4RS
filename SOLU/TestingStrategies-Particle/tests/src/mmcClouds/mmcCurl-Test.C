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
    Test the MMCcurl model

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "Cloud.H"
#include "basicMixingPopeCloud.H"
#include "basicMixingPopeParticle.H"
#include "mmcVarSet.H"
#include "psiReactionThermo.H"
#include "mmcStatusMessage.H"
#include "MMCcurl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("MMCcurl-Test","[mixing]")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();

    // Clean out old files
    Foam::rmDir("0/");

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
    #include "createFieldsForMixingTimeTest.H"

    // Create mmcVarSet for creating the cloud
    mmcVarSet Xi(entries,mesh);

    INFO("Create the MixingPopeCloud");
    basicMixingPopeCloud mixingCloud
    (
        "mixingCloud",
        mesh,
        U,
        DEff,
        rho,
        gradRhoDEff,
        Xi,
        true,
        true
    );

    INFO("Create the mixing model");
    MMCcurl<basicMixingPopeCloud> mmcCurlMixingModel
    (
        mixingCloud.subModelProperties(),
        mixingCloud,
        Xi
    );
    
    
    // ========================================================================
    //           Add two particles to the cloud 
    // ========================================================================
    
    // Delete all particles in the cloud
    mixingCloud.clear();
    
    label icell = -1;
    label itetface = -1;
    label itetpt = -1;
        
    mesh.findCellFacePt(mesh.C()[0],icell,itetface,itetpt);
    
    // Create two particles and their corresponding eulerian data 
    basicMixingPopeParticle*  pPtr = new basicMixingPopeParticle
    (
        mesh,
        Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
        icell,
        itetface,
        itetpt    
    );
    
    scalarField Y(2,0);
    Y[0] = 0.2; // This is O2
    Y[1] = 0.8; // This is CH4
    
    // Add the mixture fraction
    pPtr->Y() = Y;
    pPtr->wt() = 0.1; // Give it some random weight
    // Check that it was set correctly
    REQUIRE(pPtr->Y()[0] == 0.2);
    
    mixingCloud.addParticle(pPtr);
    
    mesh.findCellFacePt(mesh.C()[50],icell,itetface,itetpt);
    
    basicMixingPopeParticle*  qPtr = new basicMixingPopeParticle
    (
        mesh,
        Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
        icell,
        itetface,
        itetpt    
    );
    
    
    Y[0] = 0.4; // This is O2
    Y[1] = 0.6; // This is CH4
    
    qPtr->Y() = Y;
    qPtr->wt() = 0.1;
    
    REQUIRE(qPtr->Y()[0]==0.4);
    
    mixingCloud.addParticle(qPtr);
    
    // Eulerian data fields 
        
    typedef typename mixParticleModel<basicMixingPopeCloud>::eulerianFieldData 
        eulerianData;
    
    eulerianData pData;
    pData.position() = pPtr->position();
    pData.XiR() = pPtr->XiR();
    pData.D() = D[0];
    pData.Dt() = Dt[0];
    pData.DeltaE() = DeltaE[0];
    
    eulerianData qData;
    qData.position() = qPtr->position();
    qData.XiR() = pPtr->XiR();
    qData.D() = D[0];
    qData.Dt() = Dt[0];
    qData.DeltaE() = DeltaE[0];
    
    // Check that the mesh size is correct:
    INFO("Check that the LES cell size matches for the given test");
    REQUIRE(Catch::Approx(DeltaE[0]) == 0.0002);
    
    // Large mixing time to see an effect
    scalar deltaT = 1E-4;
    
    mmcCurlMixingModel.mixpair
    (
        *pPtr,
        pData,
        *qPtr,
        qData,
        deltaT        // time step
    );
    
    // ========================================================================
    //                          Check the result
    // ========================================================================
    // With the given data the result can be calculated analytically 
    // for the aISO mmcCurl model
    // The CE value is set to the default of 0.1
    // The result can easily be calculated by hand
    
    
    // Distance dx() should be equal for both particles
    REQUIRE(pPtr->dx() == qPtr->dx());
    REQUIRE(pPtr->dx() == 0.001);
    
    // Resulting mass fraction for O2 in Particle 1:
    const scalar Y_O2_p = 0.20506711;
    const scalar Y_O2_q = 0.39493289;
    
    // Check the mixing of the species mass fraction
    INFO("Mass fraction of the particle p");
    REQUIRE(Catch::Approx(pPtr->Y()[0]) == Y_O2_p);
    INFO("Mass fraction of the particle q");
    REQUIRE(Catch::Approx(qPtr->Y()[0]) == Y_O2_q);
    
    // Check that the mass is conserved
    REQUIRE(Catch::Approx(pPtr->Y()[0]+pPtr->Y()[1]) == 1.0);
    REQUIRE(Catch::Approx(qPtr->Y()[0]+qPtr->Y()[1]) == 1.0);
    
    // Write out the position of the particles for visualization
    mixingCloud.writePositions();
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

