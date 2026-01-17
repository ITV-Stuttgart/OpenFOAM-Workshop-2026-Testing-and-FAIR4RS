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
    Test the evaporation rate of a n-heptane droplet

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "globalFoamArgs.H"

#include "basicEvapDropletCloud.H"
#include "AbramzonSirignanoModel.H"
#include "psiReactionThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Generate the default volFields 
namespace nHeptaneTest
{
void writeVolFields
(
    const Foam::fvMesh& mesh
);
}

TEST_CASE("nHeptaneDroplet","[serial]")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object

    // Create the mesh
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

    // Generate the required fields for the thermo object
    nHeptaneTest::writeVolFields(mesh);

    // Construct the thermo object
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();

    volScalarField rho("rho",thermo.rho());

    // Check that it uses ideal gas
    const scalar rhoIdeal = 1E+5/(8.314/28.0*1000*555);
    REQUIRE_THAT(rho[0],Catch::Matchers::WithinRel(rhoIdeal,1E-3));

    // Package the thermo object in a SLG thermo
    SLGThermo slgThermo(mesh,thermo);

    // Load the velocity field
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    dimensionedVector g(word("g"),dimVelocity/dimTime,pTraits<vector>::zero);

    // Build the Cloud
    basicEvapDropletCloud cloud("evapCloud",rho,U,g,slgThermo);


    // Add a particle

        label icell, itetface, itetpt;
        // Find the position in the center of the domain
        mesh.findCellFacePt(mesh.bounds().centre(),icell,itetface,itetpt);
        
        // Create a particle
        basicEvapDropletParcel* dropletParticlePtr = 
        new basicEvapDropletParcel
        (
            mesh,
            Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
            icell,
            itetface,
            itetpt    
        );

        dropletParticlePtr->nParticle() = 1;
        dropletParticlePtr->T() = 300;
        dropletParticlePtr->d() = 7E-4;
        dropletParticlePtr->rho() = 684;
        dropletParticlePtr->Cp() = 1000;
        dropletParticlePtr->mass0() = dropletParticlePtr->mass();
        dropletParticlePtr->U() = vector(0,0,0);

        cloud.addParticle(dropletParticlePtr);

        // Check that everything is set correctly
        REQUIRE(cloud.first()->d() == 7E-4);
        REQUIRE(cloud.first()->T() == 300);
        REQUIRE(cloud.first()->rho() == 684);
        REQUIRE(cloud.first()->U().x() == 0);
        REQUIRE(cloud.first()->U().y() == 0);
        REQUIRE(cloud.first()->U().z() == 0);

        
        typename basicEvapDropletParcel::trackingData td(cloud);
        
        // setCellValues() is normally called in KinematicParcel::move()
        dropletParticlePtr->setCellValues(cloud,td);

    // Time Loop

        // Sets the trackTime in cloudSolution
        cloud.solution().canEvolve();

        const scalar deltaT = 1E-3;

        runTime.setDeltaT(deltaT);

        std::vector<FixedList<scalar,3>> timeEvolutionDiameter;
        timeEvolutionDiameter.reserve(5.0/deltaT);

        for (scalar t0=0; t0 < 5; t0 +=deltaT)
        {
            runTime++;
            cloud.first()->calc(cloud,td,deltaT);
            FixedList<scalar,3> temp;
            temp[0] = t0;
            temp[1] = cloud.first()->d();
            temp[2] = cloud.first()->T();
            timeEvolutionDiameter.emplace_back
            (
                std::move(temp)
            );
            if (td.keepParticle == false)
                break;
        }

        // Make a check if the evaporation rate is within
        // reasonable bounds
        // Divide the time by d0^2 to make a direct comparison
        // to the plots generated with plotDroplet.ipynb
        // Note: unit is s/mm^2 not meter squared
        REQUIRE_THAT
        (
            timeEvolutionDiameter.back()[0]/std::pow(7E-4*1000.0,2),
            Catch::Matchers::WithinRel(5,0.1)
        );

        // Write the droplet diameter to a file
        std::ofstream of("dropletDiameter_nHeptane.dat");
        for (auto& e : timeEvolutionDiameter)
            of << e[0]<<"\t"<<e[1]<<"\t"<<e[2]<<std::endl;
        of.close();
}


void nHeptaneTest::writeVolFields
(
    const Foam::fvMesh& mesh
)
{
    auto& runTime = mesh.time();

    // Nitrogen atmosphere
    volScalarField N2
    (
        IOobject
        (
            "N2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho0",dimless,1.0)
    );
    N2.write();

    // YDefault
    volScalarField Ydefault
    (
        IOobject
        (
            "Ydefault",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho0",dimless,0.0)
    );
    Ydefault.write();

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("p",dimPressure,1E+5)
    );
    p.write();

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("T",dimTemperature,555)
    );
    T.write();

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<vector>("u",dimVelocity,pTraits<vector>::zero)
    );
    U.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

