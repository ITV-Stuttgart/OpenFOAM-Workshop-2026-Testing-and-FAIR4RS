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
    Test k-d tree with construction of given points 

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "fvCFD.H"
#include "kdTree.H"
#include <random>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class VectorType>
Foam::label bruteForceSearch(const Foam::List<VectorType>&, const VectorType&);

TEST_CASE("kdTree Simple Test","[kdTree][mmcSupport]")
{
    SECTION("Only x coordinates")
    {
        // Create a list of points
        // Random number generator of C++
        std::random_device rd;
        std::mt19937 mt(rd());
        
        std::uniform_real_distribution<double> randx(0.0,10.0);
        
        const label nPoints = 10000;
        List<vector> points(nPoints);
        for (vector& v : points)
        {
            v.x() = randx(mt);
            v.y() = 0;
            v.z() = 0;
        }
        
        // Create k-d tree 
        List<scalar> weights(3,1.0);
        kdTree<vector> tree
        (
            points,
            weights,
            true,
            false
        );
        
        // Search closest points
        label nQueryPoints = 100;
        List<vector> queryPoints(nQueryPoints);
        for (vector& v : queryPoints)
        {
            v.x() = randx(mt);
            v.y() = 0;
            v.z() = 0;
        }
        
        for (auto& q : queryPoints)
        {
            auto res = tree.nNearest(q,1);
            
            REQUIRE(res[0].idx == bruteForceSearch(points,q));
        }
    }

    SECTION("x and y coordinates")
    {
        // Create a list of points
        // Random number generator of C++
        std::random_device rd;
        std::mt19937 mt(rd());
        
        std::uniform_real_distribution<double> randx(0.0,10.0);
        std::uniform_real_distribution<double> randy(0.0,1.0);
        
        const label nPoints = 10000;
        List<vector> points(nPoints);
        for (vector& v : points)
        {
            v.x() = randx(mt);
            v.y() = randy(mt);
            v.z() = 0;
        }
        
        // Create k-d tree 
        List<scalar> weights(3,1.0);
        kdTree<vector> tree
        (
            points,
            weights,
            true,
            false
        );
        
        // Search closest points
        label nQueryPoints = 100;
        List<vector> queryPoints(nQueryPoints);
        for (vector& v : queryPoints)
        {
            v.x() = randx(mt);
            v.y() = randy(mt);
            v.z() = 0;
        }
        
        for (auto& q : queryPoints)
        {
            auto res = tree.nNearest(q,1);
            
            REQUIRE(res[0].idx == bruteForceSearch(points,q));
        }
    }
    SECTION("All coordinates")
    {
        // Create a list of points
        // Random number generator of C++
        std::random_device rd;
        std::mt19937 mt(rd());
        
        std::uniform_real_distribution<double> randx(0.0,10.0);
        std::uniform_real_distribution<double> randy(0.0,1.0);
        std::uniform_real_distribution<double> randz(0.0,5.0);
        
        const label nPoints = 10000;
        List<vector> points(nPoints);
        for (vector& v : points)
        {
            v.x() = randx(mt);
            v.y() = randy(mt);
            v.z() = randz(mt);
        }
        
        // Create k-d tree 
        List<scalar> weights(3,1.0);
        kdTree<vector> tree
        (
            points,
            weights,
            true,
            false
        );
        
        // Search closest points
        label nQueryPoints = 100;
        List<vector> queryPoints(nQueryPoints);
        for (vector& v : queryPoints)
        {
            v.x() = randx(mt);
            v.y() = randy(mt);
            v.z() = randz(mt);
        }
        
        for (auto& q : queryPoints)
        {
            auto res = tree.nNearest(q,1);
            
            REQUIRE(res[0].idx == bruteForceSearch(points,q));
        }
    }
}

TEST_CASE("kdTree Simple Test 2D","[kdTree][primitives]")
{
    // Create a list of points
    // Random number generator of C++
    std::random_device rd;
    std::mt19937 mt(rd());
    
    std::uniform_real_distribution<double> rand01(0.0,10.0);
    
    const label nPoints = 10000;
    List<vector2D> points(nPoints);
    for (auto& v : points)
    {
        v.x() = rand01(mt);
        v.y() = rand01(mt);
    }
    
    // Create k-d tree 
    List<scalar> weights(2,1.0);
    kdTree<vector2D> tree
    (
        points,
        weights,
        false,
        false
    );
    
    // Search closest points
    label nQueryPoints = 100;
    List<vector2D> queryPoints(nQueryPoints);
    for (auto& v : queryPoints)
    {
        v.x() = rand01(mt);
        v.y() = rand01(mt);
    }
    
    for (auto& q : queryPoints)
    {
        auto res = tree.nNearest(q,1);
        REQUIRE(res[0].idx == bruteForceSearch(points,q));
    }
}


template<class VectorType>
Foam::label bruteForceSearch(const Foam::List<VectorType>& points, const VectorType& q)
{
    // Loop over all points
    scalar minDist = 1E+30;
    label ind = -1;
    forAll(points,i)
    {
        const auto& p = points[i];
        if (mag(p-q) < minDist)
        {
            minDist = mag(p-q);
            ind = i;
        }
        
    }
    
    return ind;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

