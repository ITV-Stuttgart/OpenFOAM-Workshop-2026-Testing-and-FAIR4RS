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

#include "particlePtrList.H"
#include <iostream>
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("particlePtrList Test","[mmcSupport]")
{
	SECTION("Constructor")
	{
		SECTION("Default Constructor")
		{
			// Default empty constructor
			Foam::particlePtrList<int> ptrList;
			REQUIRE(ptrList.size() == 0);
		}
		SECTION("Constructor with given size")
		{
			// Create with given size
			Foam::particlePtrList<int> ptrList(100);
			REQUIRE(ptrList.size() == 100);
		}
	}
	SECTION("Element Access")
	{
			// Create with given size
			Foam::particlePtrList<int> ptrList(10);
			REQUIRE(ptrList.size() == 10);
			
			// create a list of 10 integer pointers
			int* intPtr[10];
			for (int i=0; i < 10; i++ )
				intPtr[i] = new int(i);
				
			// Populate the list
			for (int i=0; i < 10; i++ )
				ptrList.set(intPtr[i],i);
			
			INFO("Loop with index");
			for (int i=0; i < 10; i++ )
				REQUIRE(ptrList[i] == i);
			
				
			// Loop over list and print result using auto range
			INFO("Loop with auto range");
			int k=0;
			for (auto e : ptrList)
				REQUIRE(e == k++);
				
			// Clear list and check if it affects the original data
			ptrList.clear();
			for (int i=0; i < 10; i++ )
				REQUIRE(*intPtr[i] == i);
			
			
			
	}		
	
	
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

