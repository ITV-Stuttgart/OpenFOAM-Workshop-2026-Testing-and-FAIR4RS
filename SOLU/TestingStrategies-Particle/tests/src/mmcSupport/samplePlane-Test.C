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
    Test the samplePlane class

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "samplePlane.H"
#include "IFstream.H"
#include <iostream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("samplePlane Test","[mmcSupport]")
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

    IOdictionary cloudProperties
    (
        IOobject
        (
            "cloudProperties",
            runTime.constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
   
    point rotAxis(0,0,1);
    point normal(1.0,1.0,0);
    point origin(0,0,0);
    

    // Create the dictionary
    dictionary dict;
    dict.getOrAdd<vector>("rotAxis",rotAxis);
    dict.getOrAdd<word>("planeType","pointAndNormal");
    dict.getOrAdd<Switch>("triangulate",false);
    auto& subDict = dict.subDictOrAdd("pointAndNormalDict");
    subDict.getOrAdd<point>("point",origin);
    subDict.getOrAdd<point>("normal",normal);

    Info << dict << endl;

    // Create the sampleMesh
    mmcSupport::samplePlane plane(mesh,dict);
    
    // Create a field with the size of the sample mesh
    Field<scalar> field(plane.size());
    forAll(field,i)
    {
        field[i] = 0;
    }
    
    // Transform point coordinates in axial and radial
    // to check with the plane
    auto transform = [&](const vector& p) -> vector2D
    {
        // axial component
        const scalar ax = (p & rotAxis);
        const point pRad = p - rotAxis*ax;
        scalar r = mag(pRad);

        // third component is always zero 
        return vector2D(ax,r);
    };

    // The domain for this test is a cylinder with:
    //  D = 1.0
    //  L = 4.0
    //  Centerpoint in origin (0 0 0)

    // Take a random point
    point p(0.1,0.1,1.0);
    
    label idx = plane.cellIndex(p);
    field[idx] =1000;

    plane.writeField("testField",field);
    Info << "===========================================================" << nl
         << "Check if samplePlane works by checking with paraview if the" << nl
         << "point located at "<<p << " is setting the plane value to 1000" << nl
         << "===========================================================" << nl
         << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

