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
    Test the particleNumberController class 

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "particleNumberController.H"
#include "Cloud.H"
#include "basicItoPopeParticle.H"
#include "barycentric.H"
#include "Random.H"

#include <random>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("particleNumberController-Test","[particleNumberControl]")
{
    // =========================================================================
    //                      Prepare Case
    // =========================================================================
    // Replace setRootCase.H for Catch2   
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


    // Global settings
    // ===============

    // Target number of particles in each super cell
    const label Npc = 20;
    const label Nhi = 22;
    const label Nlo = 15;

    
    // Create a cloud with particles
    Cloud<basicItoPopeParticle> cloud
    (
        mesh,
        "testCloud"
    );
    
    // Random number generator of C++
    std::random_device rd;
    std::mt19937 mt(rd());
    
    std::uniform_real_distribution<double> rand01(0.0,1.0);
    
    // OpenFOAM random generator
    Random rndGen;


    dictionary pControllerDict;

    // Add the required settings to the dictionary
    pControllerDict.add<scalar>("Npc",Npc);
    pControllerDict.add<scalar>("Nlo",Nlo);
    pControllerDict.add<scalar>("Nhi",Nhi);
    
    // Test constructor
    particleNumberController pController(mesh,pControllerDict);

    auto cellsInSuperCell = pController.cellsInSuperCell();

    // Loop over each super cell and add Npc particles 
    forAll(cellsInSuperCell,superCellI)
    {
        DynamicList<label> cells = cellsInSuperCell[superCellI];
        for (label i=0; i < Npc; i++)
        {
            // Select a random cell 
            label ind = (cells.size()-1)*rand01(mt);
            label celli = cells[ind];
            vector pos
            (
                mesh.C()[celli]
            );
            
            label icell = -1;
            label itetface = -1;
            label itetpt = -1;
                
            mesh.findCellFacePt(pos,icell,itetface,itetpt);
            
            // generate a random barycentric value
            
            Foam::barycentric b = barycentric01(rndGen);
            scalar maxMass = 1.0;
            if (Pstream::parRun())
                maxMass = 10.*(Pstream::myProcNo()+1);
            
            basicItoPopeParticle* pPtr = new basicItoPopeParticle
            (
                mesh,
                b,
                icell,
                itetface,
                itetpt,
                maxMass+rand01(mt),
                1.0,
                maxMass+rand01(mt), // this is wt()
                vector(1.0,1.0,1.0)
            );

            cloud.addParticle(pPtr);
        }
    }

    cloud.writePositions();
    
    // Check how many particles there are in each super cell
    for (label superCellI=0; superCellI < pController.nSuperCells(); superCellI++)
    {
        // List of particle pointers in the super cell
        auto particlePtr = pController.getParticlesInSuperCell(cloud,superCellI);
        REQUIRE(particlePtr.size() == Npc);
    }

    // =========================================================================
    //                              Start Test
    // =========================================================================

    SECTION("Add too many particles")
    {
        // Add again Npc particles
        forAll(cellsInSuperCell,superCellI)
        {
            DynamicList<label> cells = cellsInSuperCell[superCellI];
            // Add twice the difference from Nhi-Npc
            for (label i=0; i < 2.0*(Nhi-Npc); i++)
            {
                // Select a random cell 
                label ind = (cells.size()-1)*rand01(mt);
                label celli = cells[ind];
                vector pos
                (
                    mesh.C()[celli]
                );
                
                label icell = -1;
                label itetface = -1;
                label itetpt = -1;
                    
                mesh.findCellFacePt(pos,icell,itetface,itetpt);
                
                // generate a random barycentric value
                
                Foam::barycentric b = barycentric01(rndGen);
                scalar maxMass = 1.0;
                if (Pstream::parRun())
                    maxMass = 10.*(Pstream::myProcNo()+1);
                
                basicItoPopeParticle* pPtr = new basicItoPopeParticle
                (
                    mesh,
                    b,
                    icell,
                    itetface,
                    itetpt,
                    maxMass+rand01(mt),
                    1.0,
                    maxMass+rand01(mt), // this is wt()
                    vector(1.0,1.0,1.0)
                );

                cloud.addParticle(pPtr);
            }
        }

        for (label superCellI=0; superCellI < pController.nSuperCells(); superCellI++)
        {
            // List of particle pointers in the super cell
            auto particlePtr = pController.getParticlesInSuperCell(cloud,superCellI);
            INFO("Ensure that twice as many particles per super cell exist");
            REQUIRE(particlePtr.size() == Npc+2*(Nhi-Npc));
        }

        // Correct the particle cloud
        pController.correct(cloud);

        for (label superCellI=0; superCellI < pController.nSuperCells(); superCellI++)
        {
            // List of particle pointers in the super cell
            auto particlePtr = pController.getParticlesInSuperCell(cloud,superCellI);
            
            REQUIRE(particlePtr.size() == Npc);
        }
    }

    SECTION("Too few particles")
    {
        // Add again Npc particles
        forAll(cellsInSuperCell,superCellI)
        {
            auto particlePtr = pController.getParticlesInSuperCell(cloud,superCellI);
            // Add twice the difference from Nhi-Npc
            for (label i=0; i < 2.0*(Npc-Nlo); i++)
            {
                cloud.deleteParticle(*particlePtr[i]);
            }
        }

        for (label superCellI=0; superCellI < pController.nSuperCells(); superCellI++)
        {
            // List of particle pointers in the super cell
            auto particlePtr = pController.getParticlesInSuperCell(cloud,superCellI);
            INFO("Ensure that fewer than Nlo particles exist");
            REQUIRE(particlePtr.size() < Nlo);
        }

        // Correct the particle cloud
        pController.correct(cloud);

        // Ensure that all particles have unique IDs
        HashSet<label> particleIDList;

        // Loop over all particles in the cloud and check that the particle
        // ID is unique
        for (auto p : cloud)
        {
            // Require that it is not yet in the particleIDList
            INFO("Particle ID " << p.origId() << " is a duplicate");
            REQUIRE(!particleIDList.contains(p.origId()));
            particleIDList.insert(p.origId());
        }

        for (label superCellI=0; superCellI < pController.nSuperCells(); superCellI++)
        {
            // List of particle pointers in the super cell
            auto particlePtr = pController.getParticlesInSuperCell(cloud,superCellI);
            REQUIRE(particlePtr.size() == Npc);
        }
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

