/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Test the perfect gas equation of state model

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "globalFoamArgs.H"

#include "perfectGas.H"
#include "specie.H"

TEST_CASE("nHeptaneDroplet","[serial]")
{
    dictionary specieDict;
    specieDict.add("massFraction",1.0);
    specieDict.add("molWeight",28.0);

    // Add specie as sub dictionary
    dictionary dict;
    dict.add("specie",specieDict);

    perfectGas<specie> EoS(dict);

    // Check a boolean value
    REQUIRE(perfectGas<specie>::incompressible==false);

    // Check for a floating point number comparison
    REQUIRE_THAT(EoS.rho(1E+5,300),Catch::Matchers::WithinRel(1.122540,1E-6));
}