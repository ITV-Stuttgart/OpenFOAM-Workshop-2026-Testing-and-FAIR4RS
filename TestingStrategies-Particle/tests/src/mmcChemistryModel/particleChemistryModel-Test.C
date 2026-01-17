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
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

Description
    Test the particleChemistryModel class

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "psiReactionThermo.H"
#include "particleChemistryModel.H"
#include "particleTDACChemistryModel.H"
#include "thermoPhysicsTypes.H"
#include "ode.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("particleChemistryModel","[mmcChemistryModel]")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object
    #include "createMesh.H"        // create the mesh object

    Info<< "Reading thermophysical properties\n" << endl;
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");


    // Create a list of species
    const label nComponents = 7;
    scalarField Y(nComponents,0);
    Y[0] = 0;   // Species O
    Y[1] = 0.5; // Species O2
    Y[2] = 0.5; // Species H2


    using particleChemModel = ode<particleChemistryModel<psiReactionThermo,gasHThermoPhysics>>;

    particleChemModel chemModel(thermo);


    // If T < Treact no reaction should happen
    // Treact is set in the constant/chemistryProperties file 
    {
        const scalar time = chemModel.particleCalculate
        (
            0,
            1e-6,
            1,
            1e5,
            300,
            Y
        );
        REQUIRE(time==GREAT);
    }


    const scalar time = chemModel.particleCalculate
    (
        0,
        1e-6,
        1,
        1e5,
        1000,
        Y
    );

    // Value taken from old 5.x version commit: 
    // 60a6064ebca1983ced9130bccf58c64f05600d54
    REQUIRE(Catch::Approx(time)==8.28966E-6);

    REQUIRE(Catch::Approx(Y[1]) == 0.499997);
    REQUIRE(Catch::Approx(Y[2]) == 0.5);
    REQUIRE(Catch::Approx(Y[3]) == 0);
    REQUIRE(Catch::Approx(Y[4]) == 0);
    REQUIRE(Catch::Approx(Y[5]) == 3.32563E-6);
    REQUIRE(Catch::Approx(Y[6]) == 0);

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

