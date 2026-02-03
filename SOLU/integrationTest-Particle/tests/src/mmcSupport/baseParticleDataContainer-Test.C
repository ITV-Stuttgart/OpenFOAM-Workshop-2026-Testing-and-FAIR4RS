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
    Test the baseParticleDataContainer class

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "baseParticleDataContainer.H"
#include "IFstream.H"
#include <iostream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("baseParticleDataContainer Test","[mmcSupport]")
{
	SECTION("Constructor")
	{
        // Create with given size
        baseParticleDataContainer container(1,2,3);
        REQUIRE(container.position() == vector(0,0,0));

        container.position() = vector(1,2,3);
        REQUIRE(container.position() == vector(1,2,3));

        REQUIRE(container.labelValues().size() == 1);
        REQUIRE(container.scalarValues().size() == 2);
        REQUIRE(container.scalarFields().size() == 3);
        REQUIRE(container.local() == true);

        SECTION("Copy Constructor")
        {
            baseParticleDataContainer copyContainer(container);
            REQUIRE(copyContainer.position() == vector(1,2,3));
            REQUIRE(copyContainer.labelValues().size() == 1);
            REQUIRE(copyContainer.scalarValues().size() == 2);
            REQUIRE(copyContainer.scalarFields().size() == 3);
            REQUIRE(copyContainer.local() == true);

        }

        SECTION("Access Operator")
        {
            REQUIRE(container.size() == 6);
            container.labelValues()[0] = 10;
            REQUIRE(container[0] ==  1);    // Access pos[0]
            REQUIRE(container[1] ==  2);    // Access pos[1]
            REQUIRE(container[2] ==  3);    // Access pos[2]
            REQUIRE(container[3] == 10);    // Access labelValues[0]
            FatalError.throwing(true);
            REQUIRE_THROWS(container[10]);
        }
	}

    SECTION("Read and write container")
    {
        // Create a container  
        baseParticleDataContainer container(1,2,1);

        container.labelValues()[0] = 2;
        container.scalarValues()[1] = 3.3;
        
        scalarField field(10,4.2);
        container.scalarFields()[0] = field;

        REQUIRE(container.labelValues()[0] == 2);
        REQUIRE(container.scalarValues()[1] == 3.3);
        REQUIRE(container.scalarFields()[0][0] == 4.2);
        REQUIRE(container.local() == true);

        SECTION("IO Stream")
        {
            Foam::OFstream ofs("particleDataTest.dat");
            ofs << container;

            Foam::IFstream ifs("particleDataTest.dat");

            baseParticleDataContainer copyContainer(ifs);

            REQUIRE(copyContainer.labelValues().size() == container.labelValues().size());
            REQUIRE(copyContainer.scalarValues().size() == container.scalarValues().size());
            REQUIRE(copyContainer.scalarFields().size() == container.scalarFields().size());
            REQUIRE(copyContainer.local() == false);

            REQUIRE(copyContainer.labelValues()[0] == container.labelValues()[0]);
            REQUIRE(copyContainer.scalarValues()[1] == container.scalarValues()[1]);
            REQUIRE(copyContainer.scalarFields()[0][0] == container.scalarFields()[0][0]);
        }

        SECTION("Append a container")
        {
            // Create a copy
            baseParticleDataContainer copyContainer(container);

            copyContainer.append(container);
            
            REQUIRE(copyContainer.labelValues().size() == 2.0*container.labelValues().size());
            REQUIRE(copyContainer.scalarValues().size() == 2.0*container.scalarValues().size());
            REQUIRE(copyContainer.scalarFields().size() == 2.0*container.scalarFields().size());
            REQUIRE(copyContainer.local() == true);

            REQUIRE(copyContainer.labelValues()[0] == copyContainer.labelValues()[0]);
            REQUIRE(copyContainer.scalarValues()[1] == copyContainer.scalarValues()[1]);
            REQUIRE(copyContainer.scalarFields()[0][0] == copyContainer.scalarFields()[0][0]);
        }
    }


}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

