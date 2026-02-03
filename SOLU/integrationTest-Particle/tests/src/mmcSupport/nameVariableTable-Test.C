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
    Test the nameVariableTable class

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "nameVariableTable.H"

#include <map>

TEST_CASE("nameVariableTable Test","[mmcSupport]")
{
    nameVariableTable table;

    scalar rho = 10;
    scalar T   = 300;
    scalar alpha=1.0;

    std::map<Foam::word,Foam::scalar> varNameSet = 
    {
        {"rho",rho},
        {"T",T},
        {"alpha",alpha}
    };

    table.addNamedVariable("rho",rho);
    table.addNamedVariable("T",T);
    table.addNamedVariable("alpha",alpha);

    REQUIRE(table.get("rho") == rho);
    REQUIRE(table.get("T") == T);
    REQUIRE(table.get("alpha") == alpha);

    // Check that a change in the variable is reflected in the table
    T = 310;
    REQUIRE(table.get("T") == T);

    // Update out control map
    varNameSet["T"] = T;

    SECTION("Get all variables")
    {
        // Check that the names are correct:
        wordList varNames = table.getAllVarNames();
        label foundVars=0;
        for (auto& var :  varNames)
        {
            auto it = varNameSet.find(var);
            if (it != varNameSet.end())
                foundVars++;
        }
        REQUIRE(foundVars == label(varNameSet.size()));

        // Check that all variables are stored correctly
        DynamicList<scalar> container;
        table.storeAllVars(container);
        forAll(container,ind)
        {
            // Loop over all variable names, as the order is unspecified
            // in OpenFOAM HashTable
            auto it = varNameSet.find(varNames[ind]);
            REQUIRE(container[ind]==it->second);
        }
    }
    SECTION("Get vars of varList")
    {
        wordList vars = {"rho"};
        DynamicList<scalar> container;
        table.storeVarsByList(container,vars);
        REQUIRE(container[0] == rho);
    }
    

}