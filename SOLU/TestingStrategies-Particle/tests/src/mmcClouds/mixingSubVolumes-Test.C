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
    Test the mixingSubVolumes class

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "mixingSubVolumes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("mixingSubVolumes-Test","[mixing-parallel]")
{
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

    // Only run in parallel 
    if (!Pstream::parRun())
        FatalError << "Test must be run in parallel" << exit(FatalError);

    // create a dictionary 
    dictionary dict;

    dict.lookupOrAddDefault("maxProcInSubVolume",8);
    dict.lookupOrAddDefault("nSubVolumeSets",10);


    // Factors for normalization -- here some random values are chosen
    const scalar ri = 0.003;
    const scalar Xii = 0.003;

    // Test the constructor
    mixingSubVolumes mixSubVolumes(mesh,dict,ri,Xii);

    // Return list of sub-volumes
    auto subVolumeSets = mixSubVolumes.subVolumeSets();

    // Write subVolumeSets to visualize in paraview
    mixSubVolumes.writeSubVolumeSets();


    for (auto& subVolumeSet : subVolumeSets)
    {
        HashSet<label> foundProcessor;
        // Check that all processors are included in the set
        for (auto& subVolume : subVolumeSet)
        {
            for (auto procI : subVolume)
            {
                bool inserted = foundProcessor.insert(procI);
                INFO("Check that no processor is included twice in a sub-volume");
                INFO("Found: "+Foam::name(procI)+" twice");
                REQUIRE(inserted);
            }
        }

        REQUIRE(foundProcessor.size() == Pstream::nProcs());

        // Check that all processors found the same subvolumes
        List<DynamicList<List<label>>> subVolumeSetProcs(Pstream::nProcs());
        subVolumeSetProcs[Pstream::myProcNo()] = subVolumeSet;
        Pstream::gatherList(subVolumeSetProcs);
        Pstream::scatterList(subVolumeSetProcs);

        // Loop over all processors
        for (label procI=0; procI < Pstream::nProcs(); procI++)
        {
            for (label k=0; k < subVolumeSet.size(); k++)
            {
                auto& subVolume = subVolumeSet[k];
                for (label e=0; e < subVolume.size(); e++)
                    REQUIRE(subVolumeSetProcs[procI][k][e]==subVolumeSet[k][e]);
            }
        }
    }

    // Get a random sub-volume set and check that all processors do the same
    const List<List<label>>& subVolumeSet = mixSubVolumes.getSubVolumeSet();
    // Check that all processors found the same subvolumes
    List<DynamicList<List<label>>> subVolumeSetProcs(Pstream::nProcs());
    subVolumeSetProcs[Pstream::myProcNo()] = subVolumeSet;
    Pstream::gatherList(subVolumeSetProcs);
    Pstream::scatterList(subVolumeSetProcs);

    // Loop over all processors
    for (label procI=0; procI < Pstream::nProcs(); procI++)
    {
        for (label k=0; k < subVolumeSet.size(); k++)
        {
            auto& subVolume = subVolumeSet[k];
            for (label e=0; e < subVolume.size(); e++)
                REQUIRE(subVolumeSetProcs[procI][k][e]==subVolumeSet[k][e]);
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

