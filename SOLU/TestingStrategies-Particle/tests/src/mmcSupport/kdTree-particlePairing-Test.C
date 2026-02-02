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
    along with  mmcFoam.  If not, see <http://www.gnu.org/licenses/>.

Description
    Test the kdTree search algorithm in comparison to a brute force approach
    on a real data set

    Times of construction and search are reported

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "particleMatchingAlgorithm.H"
#include "fvCFD.H"
#include "IFstream.H"
#include <chrono>
#include <random>
#include <queue>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace particlePairingTest
{
    class dataContainer
    {
        private:
            vector position_;
            List<scalar> XiR_;
            scalar TF_;

        public:

            dataContainer() = default;

            vector& position() {return position_;}
            List<scalar>& XiR() {return XiR_;}
            scalar& TF() {return TF_;}

            const vector& position() const {return position_;}
            const List<scalar>& XiR() const {return XiR_;}
            const scalar& TF() const {return TF_;}
    };

};



TEST_CASE("Particle Matching Algorithm test","[mmcSupport]")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object

    // DATA SETTINGS

        // Data for dimensions and weights taken from an example Shear Layer test
        // case
        labelList dims(4);
        dims[0] = 0;
        dims[1] = 1;
        dims[2] = 2;
        dims[3] = 3;
        
        //dims[0] = 1;
        //dims[1] = 2;
        //dims[2] = 3;
        //dims[3] = 13;
        
    // LOAD INPUT DATA

        // First load the list of particles from the case
        Foam::IFstream ifs(args.path()/"particleDataSets/particleList.dat");
        
        const DynamicList<List<scalar>> particleListFull(ifs);

        // Store the data in the dataContainer
        DynamicList<particlePairingTest::dataContainer> pData;
        pData.reserve(particleListFull.size());
        forAll(particleListFull,i)
        {
            auto& p = particleListFull[i];
            particlePairingTest::dataContainer data;
            data.position() = vector(p[dims[0]],p[dims[1]],p[dims[2]]);
            data.XiR().resize(1);
            data.XiR()[0] = p[dims[3]];
            pData.append(std::move(data));
        }


    // Find pairs 

        const scalar ri = 1E-4;
        const List<scalar> Xii = {0.03};

        // Create the dictionary to read
        IStringStream is
        (
            "particleMatchingAlgorithm   kdTreeLike;\n"
        );
        IOdictionary dict
        (
            IOobject
            (
                "dict",
                runTime.constant(),
                runTime,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            is
        );

        // Create the particle matching algorithm
        auto pMatch = 
            particleMatchingAlgorithm<particlePairingTest::dataContainer>::New
            (
                ri,
                Xii,
                dict
            );

        

        // Find the pairs and time the process
        DynamicList<List<label>> pairs;
        auto t1 = std::chrono::high_resolution_clock::now();
        pMatch->findPairs(pData,pairs);
        auto t2 = std::chrono::high_resolution_clock::now();
        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
        Info << "\tFinding pairs took: " << ms_int.count() << "ms" << endl;

        // Search the number of pairs with three entries
        label nPairsWithThreeParticles = 0;
        for (auto p : pairs)
        {
            if (p.size() > 2)
                nPairsWithThreeParticles++;
        }
        Info << "n Pairs with three particles: " << nPairsWithThreeParticles << endl;
}




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

