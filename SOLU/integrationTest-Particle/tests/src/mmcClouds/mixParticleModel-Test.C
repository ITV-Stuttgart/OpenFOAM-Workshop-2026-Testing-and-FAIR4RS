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
    Test global particle pairing model

    1. Run test in single core to create a reference solution
    2. Run test in parallel and check if the same particle pairs have been 
       found

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "Cloud.H"
#include "basicMixingPopeCloud.H"
#include "basicMixingPopeParticle.H"
#include "mmcVarSet.H"
#include "psiReactionThermo.H"
#include "mmcStatusMessage.H"
#include "dummyMixModel.H"
#include "IFstream.H"
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeCloudMixingModelType(dummyMixModel, basicMixingPopeCloud);
}


// Function to check if the particle pairs match
bool containPairInList(const List<labelList>& list, const labelList& pair)
{
    // Search through list and check if pairs match
    for (auto& compPair : list)
    {
        if (pair.size() == compPair.size())
        {
            bool containsParticle = false;
            for (auto& p1 : pair)
            {
                for (auto& p2 : compPair)
                {
                    if (p1 == p2)
                    {
                        containsParticle = true;
                        break;
                    }
                }
                if (containsParticle)
                    break;
            }
            if (containsParticle)
            {
                // check that all particles match
                for (auto& p1 : pair)
                {
                    bool foundParticleID = false;
                    for (auto& p2 : compPair)
                    {
                        if (p1 == p2)
                            foundParticleID = true;
                    }
                    if (!foundParticleID)
                        return false;
                }
                return true;
            }
        }
    }
    return false;
}


TEST_CASE("mixParticleModel-Test","[mixing]")
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

    // Set the mixture fraction field to the correct value for this test
    forAll(mesh.C(),cellI)
    {
        f[cellI] = cellI*1.0/mesh.C().size();
    }
    f.write();

    // Create mmcVarSet for creating the cloud
    mmcVarSet Xi(entries,mesh);

    // ========================================================================
    //         For Single-Core Run
    // ========================================================================
    if (!Pstream::parRun())
    {
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
        dummyMixModel<basicMixingPopeCloud> mixModel
        (
            mixingCloud.subModelProperties(),
            mixingCloud,
            Xi
        );

        // To match the particle pairs between single core and parallel cases
        // the original particle ID of the single core run is stored in the
        // weight field of the particle for this test!
        // Note: origId is reset when read in parallel

        for (auto& p : mixingCloud)
        {
            p.wt() = p.origId();
        }

        // Write out the particle position and their original ID. 
        // This allows to check if the particle positions read in parallel are
        // the same 
        OFstream of("particlePos.dat");
        of << mixingCloud.size()<<endl;
        for (auto& p : mixingCloud)
        {
            of << p.wt() <<" "<<p.XiR()[0]<<" "<<p.position()<<endl;
        }

        // ====================================================================
        //           Find particle pairs
        // ====================================================================

        INFO("Find particle pairs");
        Foam::DynamicList<Foam::labelList> pairs;
        mixModel.callFindPairs(pairs);

        // Write out particle list
        OFstream os("particlePairIDs.dat");
        os << pairs;
        os << flush;

        mixingCloud.writeFields();

        Info <<"========================================================"<<nl
             << "The single core run prepares the case for the parallel" <<nl
             << "test case. Please re-run the same test case in parallel"<<nl
             <<"========================================================"<<nl
             << endl;
    }
    // ========================================================================
    //         For Parallel runs
    // ========================================================================
    else
    {
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
        dummyMixModel<basicMixingPopeCloud> mixModel
        (
            mixingCloud.subModelProperties(),
            mixingCloud,
            Xi
        );

        // ====================================================================
        //           Check the particle positions
        // ====================================================================

        // Read particle positions of single core
        IFstream ifs("particlePos.dat");
        std::map<label,Tuple2<scalar,vector>> map;

        label origID;
        scalar pwt;
        vector pos;
        label cloudSize;
        ifs >> cloudSize;
        for (int k=0; k < cloudSize; k++)
        {
            ifs >> origID;
            ifs >> pwt;
            ifs >> pos;
            Tuple2<scalar,vector> temp(pwt,pos);
            map.insert(std::pair<label,Tuple2<scalar,vector> >(origID,temp));
        }


        // Loop over own positions
        for (auto& p : mixingCloud)
        {
            auto it = map.find(p.wt());
            if (it == map.end())
                FatalError 
                    << "Particle position not found in map"
                    << exit(FatalError);
            

            INFO("Check that particle positions of parallel and single core match");
            auto temp = it->second;
            REQUIRE_THAT(p.XiR()[0],Catch::Matchers::WithinRel(temp.first(),1E-8));
            forAll(p.position(),i)
            {
                REQUIRE_THAT(p.position()[i],Catch::Matchers::WithinRel(temp.second()[i],1E-8));
            }
        }

        // ====================================================================
        //           Find particle pairs
        // ====================================================================

        INFO("Find particle pairs");
        Foam::DynamicList<Foam::labelList> pairs;
        mixModel.callFindPairs(pairs);

        OFstream os("particlePairIDs-"+name(Pstream::myProcNo())+".dat");
        os << pairs;
        os << flush;

        // Read the found pairs from the single core solution
        IFstream is("particlePairIDs.dat");
        Foam::DynamicList<Foam::labelList> singleCorePairs;
        is >> singleCorePairs;

        for (auto& pair : pairs)
        {
            INFO
            (
                "Check pair "
                + Foam::name(pair[0]) + " " +Foam::name(pair[1])
                + " on processor " + Foam::name(Pstream::myProcNo())
            );
            REQUIRE(containPairInList(singleCorePairs,pair));
        }

    }
    return;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

