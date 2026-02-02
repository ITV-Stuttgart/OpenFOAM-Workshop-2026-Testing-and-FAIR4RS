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
    Test the particleThermoMixture function.
    Requires a test case in which the thermo model can be constructed.

    This class provides easy access to the individual species properties and 
    the mixture property of a composition. The advantage is, that the 
    composition of the field does not need to be adjusted, as would be required
    for the default mixture properties of OpenFOAM.

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2024
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "baseParticleThermoMixture.H"
#include "psiReactionThermo.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiReactionThermo> createThermoReactive
(
    const Foam::fvMesh& mesh
);


Foam::autoPtr<Foam::psiThermo> createThermoMultiComponent
(
    const Foam::fvMesh& mesh
);

TEST_CASE("particleThermoMixture-Test","[reactive]")
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


    SECTION("Multicomponent")
    {
        // Create the thermo model
        autoPtr<psiThermo> thermoPtr = createThermoMultiComponent(mesh);

        thermoPtr->correct();

        // Get the current mass fraction in cell 0
        List<scalar> Y(2,1.0);

        // Create an arbritary species composition
        forAll(Y,i)
        {
            Y[i] = 1.0/Y.size();
        }

        // Create the particle thermo mixture class
        auto particleMixturePtr = baseParticleThermoMixture::New(thermoPtr());
        particleMixturePtr->update(Y);

        // Get the pressure and temperature field
        const volScalarField& p = thermoPtr->p();
        const volScalarField& T = thermoPtr->T();

        auto rhoField = thermoPtr().rho();

        REQUIRE(particleMixturePtr->rho(p[0],T[0]) == rhoField()[0]);
    }

    SECTION("Reactive")
    {
        // Create the thermo model
        auto thermoPtr = createThermoReactive(mesh);
        thermoPtr->correct();

        // Get the current mass fraction in cell 0
        PtrList<volScalarField>& YFields= thermoPtr->composition().Y();
        List<scalar> Y(YFields.size());

        forAll(YFields,i)
        {
            Y[i] = YFields[i][0];
        }

        // Check that the sum is 1.0
        REQUIRE(Foam::sum(Y) == 1.0);

        // Create an arbritary species composition
        forAll(Y,i)
        {
            Y[i] = 1.0/Y.size();
            YFields[i][0] = Y[i];
        }

        thermoPtr->correct();

        // Create the particle thermo mixture class
        auto particleMixturePtr = baseParticleThermoMixture::New(thermoPtr());
        particleMixturePtr->update(Y);

        // Get the pressure and temperature field
        const volScalarField& p = thermoPtr->p();
        const volScalarField& T = thermoPtr->T();

        auto rhoField = thermoPtr().rho();

        REQUIRE_THAT(particleMixturePtr->rho(p[0],T[0]),Catch::Matchers::WithinRel(rhoField()[0],1E-6));
    }
}



Foam::autoPtr<Foam::psiReactionThermo>
createThermoReactive(const Foam::fvMesh& mesh)
{
    {
        // Create the pressure and temperature file
        volScalarField* p = new volScalarField
        (
            IOobject
            (
                "p",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("p",dimPressure,1E+5)
        );
        p->write();
        p->store();

        volScalarField* T = new volScalarField
        (
            IOobject
            (
                "T",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("T",dimTemperature,300)
        );
        T->write();
        T->store();
    }

    // Create the thermo dictionary to read
    IStringStream is
    (
        "thermoType"
        "{"
        "    type            hePsiThermo;"
        "    mixture         reactingMixture;"
        "    transport       sutherland;"
        "    thermo          janaf;"
        "    energy          sensibleEnthalpy;"
        "    equationOfState perfectGas;"
        "    specie          specie;"
        "}"
        "\n"
        "inertSpecie N2;"
        "\n"
        "chemistryReader foamChemistryReader;"
        "foamChemistryFile \"<constant>/reactionsGRI\";"
        "foamChemistryThermoFile \"<constant>/thermo.compressibleGasGRI\";"
    );

    IOdictionary thermophysicalPropertiesDict
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    thermophysicalPropertiesDict.regIOobject::write();

    return psiReactionThermo::New(mesh);
}


Foam::autoPtr<Foam::psiThermo>
createThermoMultiComponent(const Foam::fvMesh& mesh)
{
    {
        // Create the pressure and temperature file
        volScalarField* p = new volScalarField
        (
            IOobject
            (
                "p",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("p",dimPressure,1E+5)
        );
        p->write();
        p->store();

        volScalarField* T = new volScalarField
        (
            IOobject
            (
                "T.multicomponent",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("T",dimTemperature,300)
        );
        T->write();
        T->store();

        volScalarField* H2O = new volScalarField
        (
            IOobject
            (
                "H2O.multicomponent",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("T",dimless,0.5)
        );
        H2O->write();
        H2O->store();

        volScalarField* N2 = new volScalarField
        (
            IOobject
            (
                "N2.multicomponent",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("T",dimless,0.5)
        );
        N2->write();
        N2->store();
    }

    // Create the thermo dictionary to read
    IStringStream is
    (
        "thermoType"
        "{"
        "    type            hePsiThermo;"
        "    mixture         multiComponentMixture;"
        "    transport       const;"
        "    thermo          hConst;"
        "    energy          sensibleEnthalpy;"
        "    equationOfState perfectGas;"
        "    specie          specie;"
        "}"
        "\n"
        "species         (N2 H2O);"
        "\n"
        "N2"
        "{"
        "    specie"
        "    {"
        "        molWeight   28;"
        "    }"
        "    thermodynamics"
        "    {"
        "        Cp          1000;"
        "        Hf          0;"
        "    }"
        "    transport"
        "    {"
        "        mu          1.5e-5;"
        "        Pr          0.7;"
        "    }"
        "}"
        "\n"
        "H2O"
        "{"
        "    specie"
        "    {"
        "        molWeight   18;"
        "    }"
        "    thermodynamics"
        "    {"
        "        Cp          4187;"
        "        Hf          0;"
        "    }"
        "    transport"
        "    {"
        "        mu          1e-6;"
        "        Pr          0.7;"
        "    }"
        "}"
        "inertSpecie     N2;"
        "\n"
    );

    IOdictionary thermophysicalPropertiesDict
    (
        IOobject
        (
            "thermophysicalProperties.multicomponent",
            mesh.time().constant(),
            mesh.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    thermophysicalPropertiesDict.regIOobject::write();

    return psiThermo::New(mesh,"multicomponent");
}





// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




