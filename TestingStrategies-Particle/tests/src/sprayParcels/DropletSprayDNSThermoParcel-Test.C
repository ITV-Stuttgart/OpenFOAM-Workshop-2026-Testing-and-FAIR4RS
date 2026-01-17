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
    Test the ItoPopeParticle

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "basicDropletSprayDNSThermoParcel.H"
#include "basicReactingPopeCloud.H"
#include "basicDropletSprayDNSThermoCloud.H"
#include "mmcStatusMessage.H"
#include "mmcVarSet.H"
#include "OFstream.H"
#include "IFstream.H"
#include "unitMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace DropletSprayDNSThermoParcelTest
{
    // Generates the thermophysicalPropertiesDict in constant/
    void writeThermophysicalPropertiesDict(Foam::Time& runTime);

    // Generate the cloudProperties file for the mmc cloud
    void writeCloudProperties(Foam::Time& runTime);

    // Generate the sprayCloudProperties
    void writeSprayCloudProperties_CrankNicolson(Foam::Time& runTime);

    // Generate the sprayCloudProperties
    void writeSprayCloudProperties_EulerExplicit(Foam::Time& runTime);

    // Generate the sprayCloudProperties with OpenFOAM liquid properties
    void writeSprayCloudProperties_EulerExplicit_OFliquidProperties
    (
        Foam::Time& runTime
    );

    // Generate the chemistryProperties
    void writeChemistryProperties(Foam::Time& runTime);

    // Generate the default volFields 
    void writeVolFields
    (
        const Foam::fvMesh& mesh,
        const Foam::scalar pDefault,
        const Foam::scalar TDefault,
        const Foam::vector UDefault
    );
}

TEST_CASE("DropletSprayDNSThermoParcel-Test","[sprayParcels][sprayParcels-parallel]")
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

    // Clean any existing 0 folder 
    Foam::rmDir("0/");

    DropletSprayDNSThermoParcelTest::writeThermophysicalPropertiesDict(runTime);
    DropletSprayDNSThermoParcelTest::writeCloudProperties(runTime);
    DropletSprayDNSThermoParcelTest::writeSprayCloudProperties_EulerExplicit(runTime);
    DropletSprayDNSThermoParcelTest::writeChemistryProperties(runTime);

    // END of all necessary fields for mmcCloud

    SECTION("Construtctor")
    {
        INFO("Create DropletSprayDNSThermoParcel");
        
        label icell = -1;
        label itetface = -1;
        label itetpt = -1;
            
        mesh.findCellFacePt(mesh.C()[0],icell,itetface,itetpt);
        
        // Create a particle
        basicDropletSprayDNSThermoParcel p
        (
            mesh,
            Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
            icell,
            itetface,
            itetpt    
        );
        
        // Set particle properties 
        p.T() = 300;
        p.Cp() = 1024;
        p.rho() = 790;  // Stored in KinematicParcel
        
        
        SECTION("Particle Statistics")
        {
            // Empty list -- sample all variables from initStatiscalData
            wordList vars;
            auto statData = p.getStatisticalData(vars);
            auto statNames = p.getStatisticalDataNames(vars);
            bool foundMDot = false;
            forAll(statNames,i)
            {
                if (statNames[i] == "mDot")
                {
                    foundMDot = true;
                    INFO("mDot is zero initialized");
                    REQUIRE(statData[i] == 0);
                }
            }
            REQUIRE(foundMDot);
            
        }


        SECTION("Particle IO - ASCII")
        {
            if (!Pstream::parRun())
            {
                // Write the particle to a file
                OFstream os("DropletSprayDNSThermoParcel.txt",OFstream::streamFormat::ASCII);
                os << p;
                os << flush;
                // Read from file 
                IFstream is("DropletSprayDNSThermoParcel.txt",OFstream::streamFormat::ASCII);
                basicDropletSprayDNSThermoParcel pNew(mesh, is);
        
                REQUIRE(pNew.T() == p.T());
                REQUIRE(pNew.Cp() == p.Cp());
                REQUIRE(pNew.rho() == p.rho());
            }
        }
        
        SECTION("Particle IO - PStream")
        {
            if (Pstream::parRun())
            {
                Info << "Running parallel tests"<<endl;
                // Each processor sends now its particle to the others
                PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
                
                for (label procI=0; procI < Pstream::nProcs(); procI++)
                {
                    if (procI != Pstream::myProcNo())
                    {
                        UOPstream toBuffer(procI,pBufs);
                        toBuffer << p;
                    }
                }
                
                pBufs.finishedSends();
                
                for (label procI=0; procI < Pstream::nProcs(); procI++)
                {
                    if (procI != Pstream::myProcNo())
                    {
                        UIPstream fromBuffer(procI,pBufs);
                        basicDropletSprayDNSThermoParcel pNew(mesh,fromBuffer);
                        REQUIRE(pNew.T() == p.T());
                        REQUIRE(pNew.Cp() == p.Cp());
                        REQUIRE(pNew.rho() == p.rho());
                    }
                }
                
            }
        }
    }
}


TEST_CASE("DropletSprayDNSThermoParcel-Test ethanol droplet","[sprayParcels]")
{
    // Ethanol droplet evaporation
    // Based on: 
    // P. Narasu, S. Boschmann, P. Pöschko, F. Zhao & E. Guthei
    // "Modeling and Simulation of Single Ethanol/Water
    //  Droplet Evaporation in Dry and Humid Air", 2020
    // https://doi.org/10.1080/00102202.2020.1724980
    //
    // Case with droplet diameter 0.6mm
    // Tg = 400K
    // Urel = 2m/s


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

    // Clean any existing 0 folder 
    Foam::rmDir("0/");
    
    DropletSprayDNSThermoParcelTest::writeThermophysicalPropertiesDict(runTime);
    DropletSprayDNSThermoParcelTest::writeCloudProperties(runTime);
    DropletSprayDNSThermoParcelTest::writeChemistryProperties(runTime);

    // Generate all required fields for MMC
    DropletSprayDNSThermoParcelTest::writeVolFields(mesh,1E+5,400,vector(2,0,0));

    // Modify the species for air
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
        dimensionedScalar("rho0",dimless,0.79)
    );
    N2.write();

    volScalarField O2
    (
        IOobject
        (
            "O2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho0",dimless,0.21)
    );
    O2.write();


    // Construct the thermo model
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();

    SLGThermo slgThermo(mesh,thermo);

    // Create all necessary fields for mmcCloud
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

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo.rho()
    );

    volScalarField D
    (
        IOobject
        (
            "D",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho0",dimless,1E-9)
    );

  volScalarField dpdt
  (
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("",dimPressure/dimTime,1E-6)
    );
    volScalarField Dt
    (
        IOobject
        (
            "Dt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho0",dimless,1E-9)
    );

    
    volScalarField DEff
    (
        IOobject
        (
            "DEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("rho0",dimless,1E-9)
    );


    volScalarField DeltaE
    (
        IOobject
        (
            "DeltaE",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("DeltaE", dimLength, SMALL)
    );

    if (mesh.objectRegistry::foundObject<volScalarField>("delta"))
    {
        const volScalarField& delta =
            mesh.objectRegistry::lookupObject<volScalarField>("delta");

    // DeltaE.internalField() = delta.internalField();
        DeltaE = delta;
    }


    volVectorField gradRhoDEff
    (
        IOobject
        (
            "gradRhoDEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::grad(rho * DEff)
    );


    //- Define location of file with mmcVariables definitions
    const fileName solverDir = 
                mmcStatusMessage::solverBasePath("applications/solvers/mmc/mmcFoam"); 

    const IOdictionary control
    (
        IOobject
        (
            "mmcVariablesDefinitions",
            // Location of solver
            solverDir,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary entries(control.subDict("variables"));
    mmcVarSet Xi(entries,mesh);

    basicReactingPopeCloud mmcCloud
    (
        "mmcCloud",
        mesh,
        U,
        DEff,
        rho,
        gradRhoDEff,
        Xi,
        true,
        true
    );


    dimensioned<vector> g(dimless);

    SECTION("EulerExplicit")
    {
        DropletSprayDNSThermoParcelTest::writeSprayCloudProperties_EulerExplicit(runTime);

        // Create the coupling cloud
        basicDropletSprayDNSThermoCloud sprayCloud
        (
            "sprayCloud",
            rho,
            U,
            g,
            slgThermo,
            true
        );

        sprayCloud.setMMCCloud(mmcCloud);

        // Reference value taken from NIST database
        REQUIRE_THAT(sprayCloud.liquidProperties().pSat(300),Catch::Matchers::WithinRel(8778.8,1E-6));
        // Diffusion coefficient is the experimental value from
        // E. Fueller, P. Schettler, J. Giddings,
        // "NEW METHOD FOR PREDICTION OF BINARY GAS-PHASE DIFFUSION COEFFICIENTS",
        // Note that the constant C in Eq. (4) of the paper seems off. It should
        // be 1E-2
        REQUIRE_THAT(sprayCloud.liquidProperties().Df(101325,298,28),Catch::Matchers::WithinAbs(1.35E-5,0.15E-5));

        REQUIRE_THAT(sprayCloud.liquidProperties().Lv(1E+5,320),Catch::Matchers::WithinRel(909809.67,0.1));


        // Create a particle 
        label icell = -1;
        label itetface = -1;
        label itetpt = -1;
        
        // Find the position in the center of the domain
        mesh.findCellFacePt(mesh.bounds().centre(),icell,itetface,itetpt);
        
        // Create a particle
        basicDropletSprayDNSThermoParcel* dropletParticlePtr = 
        new basicDropletSprayDNSThermoParcel
        (
            mesh,
            Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
            icell,
            itetface,
            itetpt    
        );

        dropletParticlePtr->nParticle() = 1;
        dropletParticlePtr->T() = 293.15;
        dropletParticlePtr->d() = 6E-4; // 0.7 mm
        dropletParticlePtr->rho() = sprayCloud.liquidProperties().rhoL
        (
            1E+5,dropletParticlePtr->T()
        );
        
        sprayCloud.addParticle(dropletParticlePtr);

        REQUIRE(dropletParticlePtr->d() == 6E-4);
        REQUIRE(dropletParticlePtr->T() == 293.15);
        
        typename basicDropletSprayDNSThermoParcel::trackingData td(sprayCloud);
        
        // setCellValues() is normally called in KinematicParcel::move()
        dropletParticlePtr->setCellValues(sprayCloud,td);


        const scalar deltaT = 5E-5;

        std::vector<FixedList<scalar,3>> timeEvolutionDiameter;
        timeEvolutionDiameter.reserve(20.0/deltaT);

        // Write the clouds before they are modified
        runTime.writeNow();
        sprayCloud.write();

        sprayCloud.solution().canEvolve();
        for (scalar t0=0; t0 < 20; t0 +=deltaT)
        {
            runTime++;

            dropletParticlePtr->calc(sprayCloud,td,deltaT);
            FixedList<scalar,3> temp;
            temp[0] = t0;
            temp[1] = dropletParticlePtr->d();
            temp[2] = dropletParticlePtr->T();
            timeEvolutionDiameter.emplace_back
            (
                std::move(temp)
            );
            if (td.keepParticle == false || dropletParticlePtr->d() < 7.6e-07)
                break;
        }

        // Write the droplet diameter to a file 
        std::ofstream of("dropletDiameter_DNS_Ethanol_EulerExplicit.dat");
        for (auto& e : timeEvolutionDiameter)
            of << e[0]<<"\t"<<e[1]<<"\t"<<e[2]<<std::endl;
        of.close();


        // Write out the mmcCloud
        runTime.writeNow();
        sprayCloud.write();
    }

    SECTION("EulerExplicit OFFluidProperties")
    {
        DropletSprayDNSThermoParcelTest::writeSprayCloudProperties_EulerExplicit_OFliquidProperties(runTime);

        // Create the coupling cloud
        basicDropletSprayDNSThermoCloud sprayCloud
        (
            "sprayCloud",
            rho,
            U,
            g,
            slgThermo,
            true
        );

        sprayCloud.setMMCCloud(mmcCloud);

        // Reference value taken from NIST database
        REQUIRE_THAT(sprayCloud.liquidProperties().pSat(300),Catch::Matchers::WithinRel(8778.8,1E-6));

        // OpenFOAM uses a different function to calculate the diffusion coefficient

        REQUIRE_THAT(sprayCloud.liquidProperties().Lv(1E+5,320),Catch::Matchers::WithinRel(909809.67,0.1));


        // Create a particle 
        label icell = -1;
        label itetface = -1;
        label itetpt = -1;
        
        // Find the position in the center of the domain
        mesh.findCellFacePt(mesh.bounds().centre(),icell,itetface,itetpt);
        
        // Create a particle
        basicDropletSprayDNSThermoParcel* dropletParticlePtr = 
        new basicDropletSprayDNSThermoParcel
        (
            mesh,
            Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
            icell,
            itetface,
            itetpt    
        );

        dropletParticlePtr->nParticle() = 1;
        dropletParticlePtr->T() = 293.15;
        dropletParticlePtr->d() = 6E-4; // 0.7 mm
        dropletParticlePtr->rho() = sprayCloud.liquidProperties().rhoL
        (
            1E+5,dropletParticlePtr->T()
        );
        
        sprayCloud.addParticle(dropletParticlePtr);

        
        typename basicDropletSprayDNSThermoParcel::trackingData td(sprayCloud);
        
        // setCellValues() is normally called in KinematicParcel::move()
        dropletParticlePtr->setCellValues(sprayCloud,td);


        const scalar deltaT = 5E-5;

        std::vector<FixedList<scalar,3>> timeEvolutionDiameter;
        timeEvolutionDiameter.reserve(20.0/deltaT);

        // Write the clouds before they are modified
        runTime.writeNow();
        sprayCloud.write();


        for (scalar t0=0; t0 < 20; t0 +=deltaT)
        {
            runTime++;

            dropletParticlePtr->calc(sprayCloud,td,deltaT);
            FixedList<scalar,3> temp;
            temp[0] = t0;
            temp[1] = dropletParticlePtr->d();
            temp[2] = dropletParticlePtr->T();
            timeEvolutionDiameter.emplace_back
            (
                std::move(temp)
            );
            if (td.keepParticle == false || dropletParticlePtr->d() < 1.3e-06)
                break;
        }

        // Write the droplet diameter to a file 
        std::ofstream of("dropletDiameter_DNS_Ethanol_EulerExplicit_OFLiquid.dat");
        for (auto& e : timeEvolutionDiameter)
            of << e[0]<<"\t"<<e[1]<<"\t"<<e[2]<<std::endl;
        of.close();


        // Write out the mmcCloud
        runTime.writeNow();
        sprayCloud.write();
    }

    SECTION("CrankNicolson")
    {
        DropletSprayDNSThermoParcelTest::writeSprayCloudProperties_CrankNicolson(runTime);

        // Create the coupling cloud
        basicDropletSprayDNSThermoCloud sprayCloud
        (
            "sprayCloud",
            rho,
            U,
            g,
            slgThermo,
            true
        );

        sprayCloud.setMMCCloud(mmcCloud);

        REQUIRE_THAT(sprayCloud.liquidProperties().pSat(300),Catch::Matchers::WithinRel(8778.8,1E-6));
        // Diffusion coefficient is the experimental value from
        // E. Fueller, P. Schettler, J. Giddings,
        // "NEW METHOD FOR PREDICTION OF BINARY GAS-PHASE DIFFUSION COEFFICIENTS",
        // Note that the constant C in Eq. (4) of the paper seems off. It should
        // be 1E-2
        REQUIRE_THAT(sprayCloud.liquidProperties().Df(101325,298,28),Catch::Matchers::WithinAbs(1.35E-5,0.15E-5));

        REQUIRE_THAT(sprayCloud.liquidProperties().Lv(1E+5,320),Catch::Matchers::WithinRel(909809.67,0.1));


        // Reference value taken from NIST database
        REQUIRE_THAT(sprayCloud.liquidProperties().pSat(300),Catch::Matchers::WithinRel(8778.8,1E-6));
        // Diffusion coefficient is the experimental value from
        // E. Fueller, P. Schettler, J. Giddings,
        // "NEW METHOD FOR PREDICTION OF BINARY GAS-PHASE DIFFUSION COEFFICIENTS",
        // Note that the constant C in Eq. (4) of the paper seems off. It should
        // be 1E-2
        REQUIRE_THAT(sprayCloud.liquidProperties().Df(101325,298,28),Catch::Matchers::WithinAbs(1.35E-5,0.15E-5));

        REQUIRE_THAT(sprayCloud.liquidProperties().Lv(1E+5,320),Catch::Matchers::WithinRel(909809.67,0.1));


        // Create a particle 
        label icell = -1;
        label itetface = -1;
        label itetpt = -1;
        
        // Find the position in the center of the domain
        mesh.findCellFacePt(mesh.bounds().centre(),icell,itetface,itetpt);
        
        // Create a particle
        basicDropletSprayDNSThermoParcel* dropletParticlePtr = 
        new basicDropletSprayDNSThermoParcel
        (
            mesh,
            Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
            icell,
            itetface,
            itetpt    
        );

        dropletParticlePtr->nParticle() = 1;
        dropletParticlePtr->T() = 293.15;
        dropletParticlePtr->d() = 6E-4; // 0.7 mm
        dropletParticlePtr->rho() = sprayCloud.liquidProperties().rhoL
        (
            1E+5,dropletParticlePtr->T()
        );
        
        sprayCloud.addParticle(dropletParticlePtr);

        
        typename basicDropletSprayDNSThermoParcel::trackingData td(sprayCloud);
        
        // setCellValues() is normally called in KinematicParcel::move()
        dropletParticlePtr->setCellValues(sprayCloud,td);


        const scalar deltaT = 5E-4;

        std::vector<FixedList<scalar,3>> timeEvolutionDiameter;
        timeEvolutionDiameter.reserve(20.0/deltaT);

        // Write the clouds before they are modified
        runTime.writeNow();
        sprayCloud.write();


        for (scalar t0=0; t0 < 20; t0 +=deltaT)
        {
            runTime++;

            dropletParticlePtr->calc(sprayCloud,td,deltaT);
            FixedList<scalar,3> temp;
            temp[0] = t0;
            temp[1] = dropletParticlePtr->d();
            temp[2] = dropletParticlePtr->T();
            timeEvolutionDiameter.emplace_back
            (
                std::move(temp)
            );
            if (td.keepParticle == false || dropletParticlePtr->d() < 4.2e-06)
                break;
        }

        // Write the droplet diameter to a file 
        std::ofstream of("dropletDiameter_DNS_Ethanol_CrankNicolson.dat");
        for (auto& e : timeEvolutionDiameter)
            of << e[0]<<"\t"<<e[1]<<"\t"<<e[2]<<std::endl;
        of.close();


        // Write out the mmcCloud
        runTime.writeNow();
        sprayCloud.write();
    }
}

// =============================================================================


void DropletSprayDNSThermoParcelTest::writeThermophysicalPropertiesDict(Foam::Time& runTime)
{
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
        ""
        "inertSpecie N2;"
        "\n"
        "chemistryReader      chemkinReader;"
        "\n"
        "// Single-step mechanism (5 species, 1 reaction)\n"
        "CHEMKINFile       \"$FOAM_CASE/constant/chemkin/chem.inp\";"
        "CHEMKINThermoFile \"$FOAM_CASE/constant/chemkin/therm.dat\";"
        "\n"
        "CHEMKINTransportFile \"$FOAM_CASE/constant/chemkin/transportProperties\";"
        "\n"
    );

    IOdictionary thermophysicalPropertiesDict
    (
        IOobject
        (
            "thermophysicalProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    thermophysicalPropertiesDict.regIOobject::write();
}


void DropletSprayDNSThermoParcelTest::writeCloudProperties(Foam::Time& runTime)
{
    IStringStream is
    (
        "solution\n"
        "{"
        "    interpolationSchemes"
        "    {"
        "        U           cellPoint;"
        "        DEff        cellPoint;"
        "        rho         cellPoint;"
        "        gradRho     cellPoint;"
        "        gradRhoDEff cellPoint;"
        "        f           cellPoint;"
        "        p           cellPoint;"
        "        dpdt        cellPoint;"
        "    }"
        "}"
        "\n"
        "particleManagement"
        "{"
        "    numCtlOn         on;"
        "\n"
        "    Npc 10;"
        "    Nlo 8;"
        "    Nhi 14;"
        "}"
        "\n"
        "subModels"
        "{"
        "    InflowBoundaryModel          NoInflow;"
        "\n"
        "    compositionModel singlePhaseMixture;"
        "\n"
        "    radiationModel OpticallyThinRadiation;"
        "    OpticallyThinRadiationCoeffs"
        "    {"
        "        TRef     300;"
        "    }"
        "\n"
        "    printMoleFractions false;"
        "    "
        "    thermoPhysicalCouplingModel KernelEstimation;"
        "    "
        "    KernelEstimationCoeffs"
        "    {"
        "        tauRelaxBlending    false;"
        "        tauRelaxStart       1000;"
        "        tauRelaxDelta       2;"
        "        tauRelax            10;"
        "        tauUnits            timestep;"
        "        fLow                0.03;   //lower flamability limit\n"
        "        fHigh               0.19;   //upper flamability limit\n"
        "        XiCm"
        "        {"
        "            zm              0.03;    "
        "        }"
        "        coupleSpecie (O2); //Dummy value"
        "        "
        "        condVariable        f;"
        "        dfMax               0.0075; //Maximum distance in f space allowed to do extrapolation \n"
        "        C2                  0.1;    //used to control resolution\n"
        "    }"
        "\n"
        "    sootModel       none;"
        "\n"
        "    mixingModel	    MMCcurl;"
        "\n"
        "    MMCcurlCoeffs"
        "    {"
        "        //- Parameters ri and fm as per Eq.(29) Cleary & Klimenko, Phys Fluids 23, 115102, 2011\n"
        "        ri                  6.7e-3;"
        "\n"
        "        Ximi"
        "        {"
        "            fm              0.03;    "
        "        }"
        "\n"
        "        localnessLimited    false;"
        "\n"
        "        fLowMixOnMaster     0.01;"
        "        fHighMixOnMaster    1;"
        "\n"
        "        fullSort            false;"
        "        aISO                true; "
        "        meanTimeScale       true;"
        "    }"
        "\n"
        "    reactionModel finiteRateParticleChemistry;"
        "\n"
        "    finiteRateParticleChemistryCoeffs"
        "    {"
        "        zLower              0.1;"
        "        zUpper              0.9;"
        "    }"
        "}"
        "\n"
        "\n"
        "eulerianStatistics"
        "{"
        "    enabled false;"
        "    axisymmetric false;"
        "    symmetryPlane zx;"
        "}"
        "\n"
        "particleSampling"
        "{"
        "    enabled             false;"
        "    "
        "\n"
        "    // mixture fraction range of sampled particles"
        "    zLower              0.01;"
        "    zUpper              1.0;"
        "\n"
        "    // only used in parallel mode"
        "    // false:    particle data will be gathered on master node and then be written to ascii data\n"
        "    // true:    if probe is located within one rank only, gathering will be ommited (which may gives speedup)\n"
        "    parallelSampling    false;"
        "\n"
        "    numParticlesToStore 1000;"
        "\n"
        "    particleSampleProbes"
        "    {}"
        "}"
    );

    IOdictionary cloudProperties
    (
        IOobject
        (
            "cloudProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    cloudProperties.regIOobject::write();
}


void DropletSprayDNSThermoParcelTest::writeSprayCloudProperties_EulerExplicit(Foam::Time& runTime)
{
    // The KinematicCloud of OpenFOAM requires a cloudName + "Properties" file
    // in the constant/ folder
    IStringStream is
    (
        "solution"
        "{"
        "   active          true;"
        "   coupled         yes;"
        "   cellValueSourceCorrection no;"
        "   maxCo           0.3;"
        "   calcFrequency   1;"
        "   maxTrackTime    1;"
        "\n"
        "   sourceTerms"
        "   {"
        "      schemes"
        "      {"
        "          U          semiImplicit 1;"
        "          rho        explicit 1;"
        "          Yi         explicit 1;"
        "          h          explicit 1;"
        "      }"
        "      resetOnStartup   false;"
        "   }"
        "\n"
        "   interpolationSchemes"
        "   {"
        "       rho             cell;"
        "       U               cellPoint;"
        "       thermo:mu       cell;"
        "       T               cell;"
        "       p               cell;"
        "       f               cell;"
        "       fGradSqr        cell;"
        "       Cp              cell;"
        "       kappa           cell;"
        "       C2H5OH          cell;"
        "       CO2             cell;"
        "       O2              cell;"
        "       H2O             cell;"
        "       N2              cell;" 
        "   }"
        "\n"
        "   integrationSchemes"
        "   {"
        "      U               Euler;"
        "      T               EulerExplicit;"
        "   }"
        "}"
        "\n"
        "constantProperties"
        "{"
        "   T0 315.00;"  
        "   Cp0 2460;"
        "   Tvap0 351.5;" 
        "}"
        "\n"
        "subModels"
        "{"
            "  infinityFromSP false; \n"
            "  EulerInfinityConditions true;\n"       
            "   particleForces"
            "   {"
            "       sphereDrag;"
            "   }"
            "\n"
            "   dropletToMMCModel   nearestNeighbor;"
            "\n"
            "   dispersionModel none;"
            "\n"
            "   patchInteractionModel standardWallInteraction;"
            "   stochasticCollisionModel none;"
            "   heatTransferModel        none;"
            "   surfaceFilmModel         none;"
            "   radiation                none;"
            "   LiquidPropertiesModel    UDF;"
            "\n"
            "   standardWallInteractionCoeffs"
            "   {"
            "       type            rebound;"
            "   }"
            "\n"
            "   UDFCoeffs"
            "   {"
            "       saturationPressure"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "               std::vector<scalar> a {59.796, -6595.0, -5.0474, 6.3e-07, 2.0};  // C2H5OH\n"
            "               return ( std::exp(a[0]+a[1]/x+a[2]*std::log(x)+a[3]*std::pow(x,a[4])) );\n"
            "           #};\n"
            "       }"
            "       saturationTemperature"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "                std::vector<scalar> a {3.5154e+02, 2.5342e+01, 2.2231e+00, 2.0784e-01};  // C2H5OH\n"
            "                scalar P1 = log(x/101325.0);\n"
            "                return ( a[0]+a[1]*P1+a[2]*P1*P1+a[3]*P1*P1*P1 );\n"
            "           #};\n"
            "       }"
            "       binaryDiffusionCoefficient  5.74e-05;;\n"
            "\n"
            "       latentHeat\n"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "           scalar Tc = 516.25;\n"
            "           scalar T1 = x/Tc;\n"
            "           std::vector<scalar> a {958345.091059, -0.4134, 0.75362, 0.0};  // C2H5OH\n"
            "           return ( a[0]*pow(1.0-T1,a[1]+a[2]*T1+a[3]*T1*T1) );\n"
            "           #};\n"
            "       }"
            "       liquidHeatCapacity\n"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "               scalar a_in = 3086.27;\n"
            "               scalar b_in = -13.9017;\n"
            "               scalar c_in = 0.0449525;\n"
            "               scalar d_in = -1.89534803995077e-05;\n"
            "               scalar e_in = 0;\n"
            "               scalar f_in = 0;\n"
            "               return ((((f_in*x+e_in)*x+d_in)*x +c_in)*x+b_in)*x+a_in;"
            "           #};\n"
            "       }"
            "       liquidAbsoluteEnthalpy\n"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "            std::vector<scalar> a {-6752827.25039109, 2052.57331394213, -0.060995463326749, -0.00238048724015426, 1.30130890620591e-05};  // C2H5OH\n"
            "            return ( a[0]+a[1]*x+a[2]*x*x+a[3]*x*x*x+a[4]*x*x*x*x );\n"
            "           #};\n"
            "       }"
            "       rhoL\n"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"                
            "                std::vector<scalar> a {70.1308387, 0.26395, 516.25, 0.2367, 0.0};  // C2H5OH\n"
            "                const scalar T = std::min(a[2],x);\n"
            "                return ( a[0]/pow(a[1],1.0+pow(1.0-T/a[2],a[3])) );\n"
            "           #};\n"
            "       }"
                  "       kappaL\n"
            "       {\n"
            "           type   coded;\n"
            "           code\n"            
            "           #{\n"
            "         std::vector<scalar> a {0.253, -0.000281, 0.0, 0.0, 0.0, 0.0};  // C2H5OH\n"
            "         return ( a[0]+a[1]*x+a[2]*x*x+a[3]*x*x*x+a[4]*x*x*x*x + a[5]*x*x*x*x*x);\n"
            "         #};\n"            
            "   }\n"
            "   }"
            "\n"
            "   nearestNeighborCoeffs\n"
            "   {\n"
            "       useSPProps                 false;  // false: fLES is used; true: zSP is used\n"
            "       lambda_x                     true;"
            "       lambda_f                    false;"
            "       lambda_T                     true;"
            "       ri                    2.5856e-04;"
            "       fm                        0.0072;"
            "       Tm                          40.0;"
            "   }"
            "   dropletComponent C2H5OH;\n"
            "   sourceDistribution  PSI;\n"
            "  distributedCoeffs\n" 
            "  {\n"
            "  type       box;\n"
            "  length     450E-6;\n"
            "  }\n"     
        "}"
    );

    IOdictionary sprayCloudProperties
    (
        IOobject
        (
            "sprayCloudProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    sprayCloudProperties.regIOobject::write();
}


void DropletSprayDNSThermoParcelTest::writeSprayCloudProperties_CrankNicolson(Foam::Time& runTime)
{
    // The KinematicCloud of OpenFOAM requires a cloudName + "Properties" file
    // in the constant/ folder
    IStringStream is
    (
        "solution"
        "{"
        "   active          true;"
        "   coupled         yes;"
        "   cellValueSourceCorrection no;"
        "   maxCo           0.3;"
        "   calcFrequency   1;"
        "   maxTrackTime    1;"
        "\n"
        "   sourceTerms"
        "   {"
        "      schemes"
        "      {"
        "          U          semiImplicit 1;"
        "          rho        explicit 1;"        
        "          Yi         explicit 1;"
        "          h          explicit 1;"
        "      }"
        "      resetOnStartup   false;"
        "   }"
        "\n"
        "   interpolationSchemes"
        "   {"
        "       rho             cell;"
        "       U               cellPoint;"
        "       thermo:mu       cell;"
        "       T               cell;"
        "       p               cell;"
        "       f               cell;"
        "       fGradSqr        cell;"
        "       Cp              cell;"
        "       kappa           cell;"
        "       C2H5OH          cell;"
        "       CO2             cell;"
        "       O2              cell;"
        "       H2O             cell;"
        "       N2              cell;" 
        "   }"
        "\n"
        "   integrationSchemes"
        "   {"
        "      U               Euler;"
        "      T               CrankNicolson;"
        "   }"
        "}"
        "\n"
        "constantProperties"
        "{"
        "   T0 315.00;"  
        "   Cp0 2460;"
        "   Tvap0 351.5;" 
        "}"
        "\n"
        "subModels"
        "{"
            "  infinityFromSP false; \n"
            "  EulerInfinityConditions true;\n"      
            "   particleForces"
            "   {"
            "       sphereDrag;"
            "   }"
            "\n"
            "   dropletToMMCModel   nearestNeighbor;"
            "\n"
            "   dispersionModel none;"
            "\n"
            "   patchInteractionModel standardWallInteraction;"
            "   stochasticCollisionModel none;"
            "   heatTransferModel        none;"
            "   surfaceFilmModel         none;"
            "   radiation                none;"
            "   LiquidPropertiesModel    UDF;"
            "\n"
            "   standardWallInteractionCoeffs"
            "   {"
            "       type            rebound;"
            "   }"
            "\n"
            "   UDFCoeffs"
            "   {"
            "       saturationPressure"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "               std::vector<scalar> a {59.796, -6595.0, -5.0474, 6.3e-07, 2.0};  // C2H5OH\n"
            "               return ( std::exp(a[0]+a[1]/x+a[2]*std::log(x)+a[3]*std::pow(x,a[4])) );\n"
            "           #};\n"
            "       }"
            "       saturationTemperature"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "                std::vector<scalar> a {3.5154e+02, 2.5342e+01, 2.2231e+00, 2.0784e-01};  // C2H5OH\n"
            "                scalar P1 = log(x/101325.0);\n"
            "                return ( a[0]+a[1]*P1+a[2]*P1*P1+a[3]*P1*P1*P1 );\n"
            "           #};\n"
            "       }"
            "       binaryDiffusionCoefficient  5.74e-05;\n"
            "\n"
            "       latentHeat\n"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "           scalar Tc = 516.25;\n"
            "           scalar T1 = x/Tc;\n"
            "           std::vector<scalar> a {958345.091059, -0.4134, 0.75362, 0.0};  // C2H5OH\n"
            "           return ( a[0]*pow(1.0-T1,a[1]+a[2]*T1+a[3]*T1*T1) );\n"
            "           #};\n"
            "       }"
            "       liquidHeatCapacity\n"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "               scalar a_in = 3086.27;\n"
            "               scalar b_in = -13.9017;\n"
            "               scalar c_in = 0.0449525;\n"
            "               scalar d_in = -1.89534803995077e-05;\n"
            "               scalar e_in = 0;\n"
            "               scalar f_in = 0;\n"
            "               return ((((f_in*x+e_in)*x+d_in)*x +c_in)*x+b_in)*x+a_in;"
            "           #};\n"
            "       }"
            "       liquidAbsoluteEnthalpy\n"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "            std::vector<scalar> a {-6752827.25039109, 2052.57331394213, -0.060995463326749, -0.00238048724015426, 1.30130890620591e-05};  // C2H5OH\n"
            "            return ( a[0]+a[1]*x+a[2]*x*x+a[3]*x*x*x+a[4]*x*x*x*x );\n"
            "           #};\n"
            "       }"
            "       rhoL\n"
            "       {\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"                
            "                std::vector<scalar> a {70.1308387, 0.26395, 516.25, 0.2367, 0.0};  // C2H5OH\n"
            "                const scalar T = std::min(a[2],x);\n"
            "                return ( a[0]/pow(a[1],1.0+pow(1.0-T/a[2],a[3])) );\n"
            "           #};\n"
            "       }"
            "       kappaL\n"
            "       {\n"
            "           type   coded;\n"
            "           code\n"            
            "           #{\n"
            "         std::vector<scalar> a {0.253, -0.000281, 0.0, 0.0, 0.0, 0.0};  // C2H5OH\n"
            "         return ( a[0]+a[1]*x+a[2]*x*x+a[3]*x*x*x+a[4]*x*x*x*x + a[5]*x*x*x*x*x);\n"
            "         #};\n"            
            "   }\n"
            "   }"
            "\n"
            "   nearestNeighborCoeffs\n"
            "   {\n"
            "       useSPProps                 false;  // false: fLES is used; true: zSP is used\n"
            "       lambda_x                     true;"
            "       lambda_f                    false;"
            "       lambda_T                     true;"
            "       ri                    2.5856e-04;"
            "       fm                        0.0072;"
            "       Tm                          40.0;"
            "   }"
            "   dropletComponent C2H5OH;\n"
            "   sourceDistribution  PSI;\n"
            "  distributedCoeffs\n" 
            "  {\n"
            "  type       box;\n"
            "  length     450E-6;\n"
            "  }\n"     
        "}"
    );

    IOdictionary sprayCloudProperties
    (
        IOobject
        (
            "sprayCloudProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    sprayCloudProperties.regIOobject::write();
}


void DropletSprayDNSThermoParcelTest::writeSprayCloudProperties_EulerExplicit_OFliquidProperties
(
    Foam::Time& runTime
)
{
    // The KinematicCloud of OpenFOAM requires a cloudName + "Properties" file
    // in the constant/ folder
    IStringStream is
    (
        "solution"
        "{"
        "   active          true;"
        "   coupled         yes;"
        "   cellValueSourceCorrection no;"
        "   maxCo           0.3;"
        "   calcFrequency   1;"
        "   maxTrackTime    1;"
        "\n"
        "   sourceTerms"
        "   {"
        "      schemes"
        "      {"
        "          rho        explicit 1;"
        "          U          semiImplicit 1;"
        "          Yi         explicit 1;"
        "          h          explicit 1;"
        "      }"
        "      resetOnStartup   false;"
        "   }"
        "\n"
        "   interpolationSchemes"
        "   {"
        "       rho             cell;"
        "       U               cellPoint;"
        "       thermo:mu       cell;"
        "       T               cell;"
        "       p               cell;"
        "       f               cell;"
        "       fGradSqr        cell;"
        "       Cp              cell;"
        "       kappa           cell;"
        "       C2H5OH          cell;"
        "       CO2             cell;"
        "       O2              cell;"
        "       H2O             cell;"
        "       N2              cell;" 
        "   }"
        "\n"
        "   integrationSchemes"
        "   {"
        "      U               Euler;"
        "      T               EulerExplicit;"
        "   }"
        "}"
        "\n"
        "constantProperties"
        "{"
        "   T0 315.00;"  
        "   Cp0 2460;"
        "   Tvap0 351.5;" 
        "}"
        "\n"
        "subModels"
        "{"
            "  infinityFromSP false; \n"
            "  EulerInfinityConditions true;\n"      
            "   particleForces"
            "   {"
            "       sphereDrag;"
            "   }"
            "\n"
            "   dropletToMMCModel   nearestNeighbor;"
            "\n"
            "   dispersionModel none;"
            "\n"
            "   patchInteractionModel standardWallInteraction;"
            "   stochasticCollisionModel none;"
            "   heatTransferModel        none;"
            "   surfaceFilmModel         none;"
            "   radiation                none;"
            "   LiquidPropertiesModel    OFLiquidProperties;"
            "\n"
            "   standardWallInteractionCoeffs"
            "   {"
            "       type            rebound;"
            "   }"
            "\n"
            "   OFLiquidPropertiesCoeffs"
            "   {"
            "       fluidName   C2H5OH;"
            "   }"
            "   dropletComponent C2H5OH;\n"
            "   sourceDistribution  PSI;\n"
            "\n"
            "  distributedCoeffs\n" 
            "  {\n"
            "  type       box;\n"
            "  length     450E-6;\n"
            "  }\n"     
            "   nearestNeighborCoeffs\n"
            "   {\n"
            "       useSPProps                 false;  // false: fLES is used; true: zSP is used\n"
            "       lambda_x                     true;"
            "       lambda_f                    false;"
            "       lambda_T                     true;"
            "       ri                    2.5856e-04;"
            "       fm                        0.0072;"
            "       Tm                          40.0;"
            "   }"
        "}"
    );

    IOdictionary sprayCloudProperties
    (
        IOobject
        (
            "sprayCloudProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    sprayCloudProperties.regIOobject::write();
}


void DropletSprayDNSThermoParcelTest::writeChemistryProperties(Foam::Time& runTime)
{
    IStringStream is
    (
        "chemistryType"
        "{"
        "    chemistrySolver   ode;"
        "    chemistryThermo   psi;"
        "}"
        ""
        "chemistry       on;"
        ""
        "initialChemicalTimeStep 1e-10;"
        ""
        "EulerImplicitCoeffs"
        "{"
        "    cTauChem        1;"
        "    equilibriumRateLimiter off;"
        "}"
        ""
        "odeCoeffs"
        "{"
        "    solver          rodas23;"
        "    absTol          1e-12;"
        "    relTol          1e-4;"
        "}"
        ""
        "reduction"
        "{"
        "    active          off;"
        "    method          pNone;;"
        "}"
        ""
        "tabulation"
        "{"
        "    active          off;"
        "    method          pNone;"
        "}"
    );

    IOdictionary chemistryProperties
    (
        IOobject
        (
            "chemistryProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    chemistryProperties.regIOobject::write();
}


void DropletSprayDNSThermoParcelTest::writeVolFields
(
    const Foam::fvMesh& mesh,
    const Foam::scalar pDefault,
    const Foam::scalar TDefault,
    const Foam::vector UDefault
)
{
    auto& runTime = mesh.time();
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

    volScalarField f
    (
        IOobject
        (
            "f",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("f",dimless,0.0)
    );
    f.write();

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
        dimensionedScalar("p",dimPressure,pDefault)
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
        dimensionedScalar("T",dimTemperature,TDefault)
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
        dimensioned<vector>("U",dimVelocity,UDefault)
    );
    U.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

