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
    Test the undirectedGraph class.

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "undirectedGraph.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("undirectedGraph","[uGraph][mmcSupport]")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();

    // Create the emtpy graph
    undirectedGraph graph;

    // Add a vertex with connections to two other vertices
    List<label> edges = {1,2,3};
    List<scalar> weights = {0.0,1.0,2.0};

    graph.addVertex(1,edges,weights);

    REQUIRE(graph.size() == 3);

    // Check the connections
    {
        auto& connectedVertices1 = graph.getNeighbor(1);
        REQUIRE(connectedVertices1.size() == 2);

        auto& connectedVertices2 = graph.getNeighbor(2);
        REQUIRE(connectedVertices2.size() == 1);
        REQUIRE(connectedVertices2[0].first == 1);
        REQUIRE(connectedVertices2[0].second == 1.0);

        auto& connectedVertices3 = graph.getNeighbor(3);
        REQUIRE(connectedVertices3.size() == 1);
        REQUIRE(connectedVertices3[0].first == 1);
        REQUIRE(connectedVertices3[0].second == 2.0);

        auto edges = graph.getEdges();
        REQUIRE(edges.size() == 2);
        REQUIRE(edges[0].w == 1.0);
        REQUIRE(edges[1].w == 2.0);

        REQUIRE(edges[0].p == 1);
        REQUIRE(edges[0].q == 2);
        
        REQUIRE(edges[1].p == 1);
        REQUIRE(edges[1].q == 3);
    }

    // Check that vertices are not added double
    List<label> newEdge = {1,3};
    List<scalar> newWeights = {1.0,3.0};

    graph.addVertex(2,newEdge,newWeights);

    REQUIRE(graph.size() == 3);
    // Check the connections
    {
        auto& connectedVertices1 = graph.getNeighbor(1);
        REQUIRE(connectedVertices1.size() == 2);

        auto& connectedVertices2 = graph.getNeighbor(2);
        REQUIRE(connectedVertices2.size() == 2);
        REQUIRE(connectedVertices2[0].first == 1);
        REQUIRE(connectedVertices2[0].second == 1.0);

        auto& connectedVertices3 = graph.getNeighbor(3);
        REQUIRE(connectedVertices3.size() == 2);
        REQUIRE(connectedVertices3[0].first == 1);
        REQUIRE(connectedVertices3[0].second == 2.0);

        auto edges = graph.getEdges();
        REQUIRE(edges.size() == 3);
        REQUIRE(edges[0].w == 1.0);
        REQUIRE(edges[1].w == 2.0);
        REQUIRE(edges[2].w == 3.0);
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

