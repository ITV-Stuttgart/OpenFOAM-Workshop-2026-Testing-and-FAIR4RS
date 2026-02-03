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
    Test the particlePtrList class

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "fvCFD.H"
#include "unitMesh.H"

#include "particleProbeCylinder.H"



TEST_CASE("particleProbeCylinder Test","[mmcSupport]")
{
    // Create a mesh object
    unitMesh uMesh(10);
    const fvMesh& mesh = uMesh.mesh();
    Info << "done"<<endl;
    // Create the dictionary for the particleProbe
    IStringStream is
    (
        "p1     (0.5 0.5 0);"
        "p2     (0.5 0.5 0.5);"
        "radius 0.3;"
    );

    dictionary dict(is);

    particleProbeCylinder probe(dict);



    {
        // Particle position
        const vector pos(0.3,0.3,0.2);

        // Create a particle 
        particle p(mesh,pos);
        INFO("Check that particle is inside (positive)");
        REQUIRE(probe.containsParticle(p) == true);
    }
    {
    // Particle position
    const vector pos(0.6,0.6,0.2);

    // Create a particle 
    particle p(mesh,pos);
    INFO("Check that particle is inside (negative)");
    REQUIRE(probe.containsParticle(p) == true);
    }
    {
    // Particle position
    const vector pos(0.1,0.1,0.2);

    // Create a particle 
    particle p(mesh,pos);
    INFO("Check that particle is outside radial direction");
    REQUIRE(probe.containsParticle(p) == false);
    }
    {
    // Particle position
    const vector pos(0.3,0.3,0.6);

    // Create a particle 
    particle p(mesh,pos);
    INFO("Check that particle is outside axial direction");
    REQUIRE(probe.containsParticle(p) == false);
    }
}
