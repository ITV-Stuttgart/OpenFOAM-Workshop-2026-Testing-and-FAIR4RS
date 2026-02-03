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

#include "kdTree.H"
#include "fvCFD.H"
#include "IFstream.H"
#include <chrono>
#include <random>
#include <queue>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Find nearest neighbors by brute force and returns their indices
namespace kdTreeTest
{

struct compareByDistance
{
    constexpr bool operator()
    (
        std::pair<Foam::scalar,Foam::label> const &a, 
        std::pair<Foam::scalar,Foam::label> const &b
    ) const noexcept
    {
        return a.first < b.first;
    }
};


Foam::labelList nNearestSearchBruteForce
(
    const List<scalarList>& particleList,
    const scalarList& query,
    const scalarList& weights,
    const scalarList& norm,
    const label n
);


void checkNeighbors
(
    const List<label>& neigh,
    const List<label>& neighKdTree
);

}

TEST_CASE("k-d Tree Test","[kdTree][mmcSupport]")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();

    // DATA SETTINGS

        // Data for dimensions and weights taken from an example Shear Layer test
        // case
        labelList dims(4);
        dims[0] = 1;
        dims[1] = 2;
        dims[2] = 3;
        dims[3] = 13;

        scalarList wts(4);
        wts[0] = 1;
        wts[1] = 1;
        wts[2] = 1;
        wts[3] = 0.03; // 0.03;
        
    // LOAD INPUT DATA

        // First load the list of particles from the case
        Foam::IFstream ifs(args.path()/"particleDataSets/particleList.dat");
        const DynamicList<List<scalar>> particleListFull(ifs);
        
        // Search 20 nearest neighbours for each cell 
        // Read the query points from file
        Foam::IFstream ifsQueryPoints(args.path()/"particleDataSets/queryPointList.dat");
        const List<scalarList> queryPoints(ifsQueryPoints);
        const label nQueryPoints = queryPoints.size();
        
    // Create a particle List based on the given dims
        List<List<scalar>> particleList(particleListFull.size());
        forAll(particleList,i)
        {
            List<scalar> temp(4);
            for (label k=0; k < 4; k++)
                temp[k] = particleListFull[i][dims[k]];
            particleList[i] = temp;
        }
        
        
    SECTION("Search nearest neighbors")
    {
        // CALCULATE NORMALIZATIONS

            // First calculate the min and max in each dimension and 
            // weigh it for the normalization factor
            scalarList norm(wts.size());

            scalarList minVal(wts.size(),1E+15);
            scalarList maxVal(wts.size(),-1E+15);

            forAll(particleList,i)
            {
                for (int j=0; j < wts.size(); ++j)
                {
                    minVal[j] = std::min(particleList[i][j],minVal[j]);
                    maxVal[j] = std::max(particleList[i][j],maxVal[j]);
                }

            }

            for (int i=0; i < wts.size(); ++i)
                norm[i] = wts[i]*(maxVal[i]-minVal[i]);
             
        // DATA TO STORE RESULTS
            
            // Store found indices of the search algorithms
            List<List<label>> foundNeighbourList_kdTree(nQueryPoints);
            List<List<label>> foundNeighbourList_kdTreeMedianBased(nQueryPoints);
            List<List<label>> foundNeighbourList_BruteForce(nQueryPoints);
        
        // GENERATE KD TREE SOLUTION - Average based sorting

            // Construct the kdTree
            auto t1 = std::chrono::high_resolution_clock::now();
            kdTree<List<scalar>> tree(particleList,wts);
            auto t2 = std::chrono::high_resolution_clock::now();
            
            auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            Info << "Construction of kdTree took: " << ms_int.count() << "ms" << endl;

            // List of results of the kdTree
            std::vector<kdTree<List<scalar>>::resList> kdTreeResultLists;
            kdTreeResultLists.reserve(nQueryPoints);


            t1 = std::chrono::high_resolution_clock::now();
            forAll(queryPoints,i)
            {
                auto res = tree.nNearest(queryPoints[i],20);
                kdTreeResultLists.push_back(std::move(res));
            }
            t2 = std::chrono::high_resolution_clock::now();
            auto timeSearched = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            Info << "Search nearest neighbor - kdTree: " << timeSearched.count() << "ms" << endl;

            // Gather statistics 
            label i=0;
            for (auto& resList : kdTreeResultLists)
            {
                labelList neighbourIndices(resList.size());
                label k = 0;
                // Gather statistics 
                for (auto& res : resList)
                    neighbourIndices[k++] = res.idx;
                foundNeighbourList_kdTree[i] = neighbourIndices;
                i++;
            }
            
        // GENERATE KD TREE SOLUTION - Median based sorting

            // Construct the kdTree
            t1 = std::chrono::high_resolution_clock::now();
            kdTree<List<scalar>> treeMedianBased(particleList,wts,true);
            t2 = std::chrono::high_resolution_clock::now();
            
            ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            Info << "Construction of kdTree (medianBased) took: " << ms_int.count() << "ms" << endl;

            // List of results of the kdTree
            std::vector<kdTree<List<scalar>>::resList> kdTreeResultLists2;
            kdTreeResultLists2.reserve(nQueryPoints);


            t1 = std::chrono::high_resolution_clock::now();
            forAll(queryPoints,i)
            {
                auto res = treeMedianBased.nNearest(queryPoints[i],20);
                kdTreeResultLists2.push_back(std::move(res));
            }
            t2 = std::chrono::high_resolution_clock::now();
            timeSearched = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            Info << "Search nearest neighbor - kdTree (medianBased): " << timeSearched.count() << "ms" << endl;

            // Gather statistics 
            i=0;
            for (auto& resList : kdTreeResultLists2)
            {
                labelList neighbourIndices(resList.size());
                label k = 0;
                // Gather statistics 
                for (auto& res : resList)
                    neighbourIndices[k++] = res.idx;
                foundNeighbourList_kdTreeMedianBased[i] = neighbourIndices;
                i++;
            }
            

        // BRUTE FORCE ALGORITHM

            t1 = std::chrono::high_resolution_clock::now();
            forAll(queryPoints,i)
            {
                auto neighbourIndices = 
                    kdTreeTest::nNearestSearchBruteForce
                    (
                        particleList,
                        queryPoints[i],
                        wts,
                        norm,
                        20
                    );
                foundNeighbourList_BruteForce[i] = std::move(neighbourIndices);
            }
            t2 = std::chrono::high_resolution_clock::now();
            auto timeSearchedBruteForce = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
            Info << "Time to search nNearest neighbours (brute force): " << timeSearchedBruteForce.count() << "ms" << endl;


        // CHECK THE IDS OF BRUTE FORCE AND KDTREE
        forAll(foundNeighbourList_kdTree,i)
        {
            kdTreeTest::checkNeighbors
            (
                foundNeighbourList_BruteForce[i],
                foundNeighbourList_kdTree[i]
            );
        }
        
        forAll(foundNeighbourList_kdTreeMedianBased,i)
        {
            kdTreeTest::checkNeighbors
            (
                foundNeighbourList_BruteForce[i],
                foundNeighbourList_kdTreeMedianBased[i]
            );
        }
    }
}


void kdTreeTest::checkNeighbors
(
    const List<label>& neigh,
    const List<label>& neighKdTree
)
{
    forAll(neigh,i)
    {
        REQUIRE(neigh[i] == neighKdTree[i]);
    }
}


Foam::labelList kdTreeTest::nNearestSearchBruteForce
(
    const List<scalarList>& particleList,
    const scalarList& query,
    const scalarList& weights,
    const scalarList& norm,
    const label n
)
{
    std::priority_queue
    <
        std::pair<scalar,label>,
        std::vector<std::pair<scalar,label>>,
        compareByDistance
    > nearestNeighbours;
    

    // Start looping over all entries
    forAll(particleList,i)
    {
        auto& particle = particleList[i];

        // Calculate distance to query point
        scalar dist = 0;
        for (int k = 0; k < weights.size(); ++k)
        {
            dist += Foam::sqr((query[k] - particle[k])/norm[k]);
        }
        
       
        if (nearestNeighbours.size() < n)
        {
            nearestNeighbours.emplace(dist,i);
            continue;
        }

        if (dist < nearestNeighbours.top().first)
            nearestNeighbours.emplace(dist,i);

        if (nearestNeighbours.size() > n)
            nearestNeighbours.pop(); 
    }
  
    if (nearestNeighbours.size() != n)
        WarningInFunction
            << "nearestNeighbours.size() "<<nearestNeighbours.size()<<" not equal "<<n<<endl;


    labelList neighbors(n);
    forAll(neighbors,i)
    {
        neighbors[i] = std::move(nearestNeighbours.top().second);
        nearestNeighbours.pop();
    }
    
    return neighbors; 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

