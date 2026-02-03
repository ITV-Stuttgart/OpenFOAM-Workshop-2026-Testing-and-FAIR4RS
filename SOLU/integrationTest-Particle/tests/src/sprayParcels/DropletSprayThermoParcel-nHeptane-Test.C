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
    Test the DropletSprayThermoParcel for an evaporating nHeptane droplet


Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "basicDropletSprayThermoParcel.H"
#include "basicReactingPopeCloud.H"
#include "basicDropletSprayThermoCloud.H"
#include "mmcStatusMessage.H"
#include "mmcVarSet.H"
#include "OFstream.H"
#include "IFstream.H"
#include "unitMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace nHeptaneTest
{
    // Generates the thermophysicalPropertiesDict in constant/
    void writeThermophysicalPropertiesDict(Foam::Time& runTime);

    // Generate the cloudProperties file for the mmc cloud
    void writeCloudProperties(Foam::Time& runTime);

    // Generate the sprayCloudProperties
    void writeSprayCloudProperties_EulerExplicit(Foam::Time& runTime);

    // Generate the sprayCloudProperties
    void writeSprayCloudProperties_EulerExplicit_OFLiquid(Foam::Time& runTime);

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


// Important Note:
// The nHeptane and Ethanol test case cannot be run together with the UDF 
// liquid property function. This is due to the way how the UDF function will
// lookup the correct function. 
// It will always lookup the first function with the matching name, which is
// then the saturation pressure function of the ethanol case. 
// Only when the program is completely terminated and restarted, it will look-up
// the correct nHeptane function. Even if the dictionary and entries are written
// correctly!

TEST_CASE("DropletSprayThermoParcel-Test nHeptane droplet")
{
    // Droplet evaporation compared to the test case of 
    // Nomura et al. (1996) and validated with the numerical results of 
    // Yang and Wong (2001)
    // Ambient conditions: 
    //  - pure nitrogen atmosphere
    //  - ambient temperatuer at 555K
    //  - queiscent air -- no convection
    // Droplet conditions
    //  - d = 0.7mm
    //  - T = 300K
    //  - Fluid: n-heptane
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

    nHeptaneTest::writeThermophysicalPropertiesDict(runTime);
    nHeptaneTest::writeCloudProperties(runTime);
    nHeptaneTest::writeChemistryProperties(runTime);
    // Generate all required fields for MMC
    nHeptaneTest::writeVolFields(mesh,1E+5,555,vector(0,0,0));
    // Modify the species for pure nitrogen atmosphere
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
        dimensionedScalar("rho0",dimless,0.0)
    );
    O2.write();
    // Construct the thermo model -- this one will be used in the mmcCloud
    autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
    psiReactionThermo& thermo = pThermo();
    auto& mixture = thermo.composition();
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

    // ===> End of generating fields for MMC <===
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

    // =====================================================================
    // Check if the stochastic particles have the correct properties
    for (auto& p : mmcCloud)
    {
        // Temperature check
        REQUIRE(p.T() == 555);
        label jN2 = mixture.species()["N2"];
        // Species composition check
        for (int i=0; i < mixture.Y().size(); i++)
        {
            if (i != jN2)
                REQUIRE(p.Y()[i] == 0);
            else
                REQUIRE(p.Y()[i] == 1.0);
        }
    }

    dimensioned<vector> g(dimless);

    SECTION("UDF")
    {

        // Time stepping settings
        const scalar deltaT = 1E-3;

        runTime.setDeltaT(deltaT);

        nHeptaneTest::writeSprayCloudProperties_EulerExplicit(runTime);
        // Create the coupling cloud
        basicDropletSprayThermoCloud sprayCloud
        (
            "sprayCloud_nHeptane",
            rho,
            U,
            g,
            mmcCloud.slgThermo(),
            true
        );
        sprayCloud.setMMCCloud(mmcCloud);

        // =========================================================================
        // Check the thermo models of the cloud
        REQUIRE_THAT
        (
            sprayCloud.liquidProperties().pSat(360.0),
            Catch::Matchers::WithinRel(71193.7,0.01)
        );

        REQUIRE_THAT
        (
            sprayCloud.liquidProperties().Lv(1E+5,320.0),
            Catch::Matchers::WithinRel(351179,0.01)
        );

        REQUIRE_THAT
        (
            sprayCloud.liquidProperties().Df(101325,298,28.0),
            Catch::Matchers::WithinRel(7.3352E-7,0.05)
        );

        // Check that the transient option is activated
        REQUIRE(sprayCloud.solution().transient() == true);

        // =========================================================================

        // Create a particle 
        label icell = -1;
        label itetface = -1;
        label itetpt = -1;
        
        // Find the position in the center of the domain
        mesh.findCellFacePt(mesh.bounds().centre(),icell,itetface,itetpt);
        
        // Create a particle
        basicDropletSprayThermoParcel* dropletParticlePtr =
            new basicDropletSprayThermoParcel
            (
                mesh,
                Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
                icell,
                itetface,
                itetpt    
            );

        dropletParticlePtr->nParticle() = 1;
        dropletParticlePtr->T() = 300;
        dropletParticlePtr->d() = 7E-4; // 0.7 mm
        dropletParticlePtr->rho() = sprayCloud.liquidProperties().rhoL
        (
            1E+5,dropletParticlePtr->T()
        );
        dropletParticlePtr->U() = vector(0,0,0);

        // Important: After commit ff4a562783bebfb1b5678cfe38782212774f243f
        dropletParticlePtr->mass0() = dropletParticlePtr->mass();
        
        // Add the particle to the cloud
        sprayCloud.addParticle(dropletParticlePtr);

        // Write the clouds before they are modified
        runTime.writeNow();
        mmcCloud.write();
        sprayCloud.write();

        // Check that the parcel in the cloud has the required properties
        REQUIRE(sprayCloud.first()->d() == 7E-4);
        REQUIRE(sprayCloud.first()->T() == 300);
        REQUIRE(sprayCloud.first()->rho() == 684);
        REQUIRE(sprayCloud.first()->U().x() == 0);
        REQUIRE(sprayCloud.first()->U().y() == 0);
        REQUIRE(sprayCloud.first()->U().z() == 0);
        REQUIRE(sprayCloud.first()->mass0() == sprayCloud.first()->mass());

        std::vector<FixedList<scalar,3>> timeEvolutionDiameter;
        timeEvolutionDiameter.reserve(5.0/deltaT);

        // Here we do not use cloud.evolve() as this is a quite slow operation
        // And the focus is on testing the evaporation rate
        typename basicDropletSprayThermoParcel::trackingData td(sprayCloud);
        // setCellValues() is normally called in KinematicParcel::move()
        dropletParticlePtr->setCellValues(sprayCloud,td);

        // Sets the trackTime in cloudSolution
        sprayCloud.solution().canEvolve();

        for (scalar t0=0; t0 < 5; t0 +=deltaT)
        {
            runTime++;
            sprayCloud.first()->calc(sprayCloud,td,deltaT);
            FixedList<scalar,3> temp;
            temp[0] = t0;
            temp[1] = sprayCloud.first()->d();
            temp[2] = sprayCloud.first()->T();
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

        // Write out the mmcCloud
        runTime.writeNow();
        mmcCloud.write();
        sprayCloud.write();
    }

    SECTION("OpenFOAM Properties")
    {

        // Time stepping settings
        const scalar deltaT = 1E-3;

        runTime.setDeltaT(deltaT);

        nHeptaneTest::writeSprayCloudProperties_EulerExplicit_OFLiquid(runTime);
        // Create the coupling cloud
        basicDropletSprayThermoCloud sprayCloud
        (
            "sprayCloud_nHeptane",
            rho,
            U,
            g,
            mmcCloud.slgThermo(),
            true
        );
        sprayCloud.setMMCCloud(mmcCloud);

        // =========================================================================
        // Check the thermo models of the cloud
        REQUIRE_THAT
        (
            sprayCloud.liquidProperties().pSat(360.0),
            Catch::Matchers::WithinRel(71193.7,0.01)
        );

        REQUIRE_THAT
        (
            sprayCloud.liquidProperties().Lv(1E+5,320.0),
            Catch::Matchers::WithinRel(351179,0.01)
        );

        // Note: OpenFOAM uses a different function to calculate the Df of 
        //       the species. See  APIdiffCoefFunc.H
        REQUIRE_THAT
        (
            sprayCloud.liquidProperties().Df(101325,298,28.0),
            Catch::Matchers::WithinRel(7.1E-6,0.05)
        );

        // Check that the transient option is activated
        REQUIRE(sprayCloud.solution().transient() == true);

        // =========================================================================

        // Create a particle 
        label icell = -1;
        label itetface = -1;
        label itetpt = -1;
        
        // Find the position in the center of the domain
        mesh.findCellFacePt(mesh.bounds().centre(),icell,itetface,itetpt);
        
        // Create a particle
        basicDropletSprayThermoParcel* dropletParticlePtr =
            new basicDropletSprayThermoParcel
            (
                mesh,
                Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
                icell,
                itetface,
                itetpt    
            );

        dropletParticlePtr->nParticle() = 1;
        dropletParticlePtr->T() = 300;
        dropletParticlePtr->d() = 7E-4; // 0.7 mm
        dropletParticlePtr->rho() = sprayCloud.liquidProperties().rhoL
        (
            1E+5,dropletParticlePtr->T()
        );
        dropletParticlePtr->U() = vector(0,0,0);

        // Important: After commit ff4a562783bebfb1b5678cfe38782212774f243f
        dropletParticlePtr->mass0() = dropletParticlePtr->mass();
        
        // Add the particle to the cloud
        sprayCloud.addParticle(dropletParticlePtr);

        // Write the clouds before they are modified
        runTime.writeNow();
        mmcCloud.write();
        sprayCloud.write();

        // Check that the parcel in the cloud has the required properties
        REQUIRE(sprayCloud.first()->d() == 7E-4);
        REQUIRE(sprayCloud.first()->T() == 300);
        REQUIRE(sprayCloud.first()->U().x() == 0);
        REQUIRE(sprayCloud.first()->U().y() == 0);
        REQUIRE(sprayCloud.first()->U().z() == 0);
        REQUIRE(sprayCloud.first()->mass0() == sprayCloud.first()->mass());

        std::vector<FixedList<scalar,3>> timeEvolutionDiameter;
        timeEvolutionDiameter.reserve(5.0/deltaT);

        // Here we do not use cloud.evolve() as this is a quite slow operation
        // And the focus is on testing the evaporation rate
        typename basicDropletSprayThermoParcel::trackingData td(sprayCloud);
        // setCellValues() is normally called in KinematicParcel::move()
        dropletParticlePtr->setCellValues(sprayCloud,td);

        // Sets the trackTime in cloudSolution
        sprayCloud.solution().canEvolve();

        for (scalar t0=0; t0 < 5; t0 +=deltaT)
        {
            runTime++;
            sprayCloud.first()->calc(sprayCloud,td,deltaT);
            FixedList<scalar,3> temp;
            temp[0] = t0;
            temp[1] = sprayCloud.first()->d();
            temp[2] = sprayCloud.first()->T();
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
        std::ofstream of("dropletDiameter_nHeptane_OFLiquid.dat");
        for (auto& e : timeEvolutionDiameter)
            of << e[0]<<"\t"<<e[1]<<"\t"<<e[2]<<std::endl;
        of.close();

        // Write out the mmcCloud
        runTime.writeNow();
        mmcCloud.write();
        sprayCloud.write();
    }
    
}

// =============================================================================

void nHeptaneTest::writeThermophysicalPropertiesDict(Foam::Time& runTime)
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
        "CHEMKINFile       \"$FOAM_CASE/constant/chemkin_nHeptane/chem.inp\";"
        "CHEMKINThermoFile \"$FOAM_CASE/constant/chemkin_nHeptane/therm.dat\";"
        "\n"
        "CHEMKINTransportFile \"$FOAM_CASE/constant/chemkin_nHeptane/transportProperties\";"
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


void nHeptaneTest::writeCloudProperties(Foam::Time& runTime)
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
        "        dpdt           cellPoint;"
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
        "        "
        "        coupleSpecie (O2); //Dummy value"
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


void nHeptaneTest::writeSprayCloudProperties_EulerExplicit(Foam::Time& runTime)
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
        "   transient       true;"
        "\n"
        "   sourceTerms"
        "   {"
        "      schemes"
        "      {"
        "          U               semiImplicit 1;"
        "      }"
        "      resetOnStartup   false;"
        "   }"
        "\n"
        "   interpolationSchemes"
        "   {"
        "      rho             cell;"
        "      U               cellPoint;"
        "      p               cell;"
        "      thermo:mu       cell;"
        "      T               cell;"
        "      Cp              cell;"
        "      kappa           cell;"
        "   }"
        "\n"
        "   integrationSchemes"
        "   {"
        "      U               Euler;"
        "      T               EulerExplicit;"
        "   }"
        "}"
        "\n"
        "subModels"
        "{"
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
            "           // From: Williamham, Taylor, et al., 1945	\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "               std::vector<scalar> a {4.02832, 1268.636,   -56.199};  // NC7H16\n"
            "               return ( std::pow(10,a[0]-a[1]/(x+a[2]))*1E+5);\n"
            "           #};\n"
            "       }"
            "       saturationTemperature"
            "       {"
            "           // From: Williamham, Taylor, et al., 1945	\n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "               std::vector<scalar> a {4.02832, 1268.636,   -56.199};  // NC7H16\n"
            "               return a[1]/(a[0]-std::log10(x/1E+5))-a[2];\n"
            "           #};\n"
            "       }"
            "\n"
            "       binaryDiffusionCoefficient  3.417931E-6;\n"
            "\n"
            "       latentHeat\n"
            "       {\n"
            "           // Taken from https://webbook.nist.gov/cgi/cbook.cgi?ID=C142825&Mask=4 \n"
            "           type    coded;\n"
            "           code\n"
            "           #{\n"
            "           scalar Tc = 540;\n"
            "           scalar Tr = x/Tc;\n"
            "           std::vector<scalar> a {5.3551E+5, 0.2831};  // NC7H16\n"
            "           return ( a[0]*std::exp(-a[1]*Tr)*std::pow(1-Tr,a[1]));\n"
            "           #};\n"
            "       }"
            "       liquidHeatCapacity\n"
            "       {\n"
            "           type    table;\n"
            "           values\n"
            "           (\n"
            "               (357 1907)\n"
            "               (373 1983)\n"
            "               (400 2103)\n"
            "               (434.35 2248)\n"
            "           );"
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
            "           type    constant;\n"
            "           value   684;\n"
            "       }"
            "   }"
            "\n"
            "   nearestNeighborCoeffs\n"
            "   {\n"
            "       useSPProps                 false;  // false: fLES is used; true: zSP is used\n"
            "       lambda_x                     1.0;"
            "       lambda_f                     0.0;"
            "       lambda_T                     1.0;"
            "       ri                    2.5856e-04;"
            "       fm                        0.0072;"
            "       Tm                          40.0;"
            "   }"
            "   dropletComponent NC7H16;\n"
        "}"
    );

    IOdictionary sprayCloudProperties
    (
        IOobject
        (
            "sprayCloud_nHeptaneProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    sprayCloudProperties.regIOobject::write();
}


void nHeptaneTest::writeSprayCloudProperties_EulerExplicit_OFLiquid(Foam::Time& runTime)
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
        "   transient       true;"
        "\n"
        "   sourceTerms"
        "   {"
        "      schemes"
        "      {"
        "          U               semiImplicit 1;"
        "      }"
        "      resetOnStartup   false;"
        "   }"
        "\n"
        "   interpolationSchemes"
        "   {"
        "      rho             cell;"
        "      U               cellPoint;"
        "      p               cell;"
        "      thermo:mu       cell;"
        "      T               cell;"
        "      Cp              cell;"
        "      kappa           cell;"
        "   }"
        "\n"
        "   integrationSchemes"
        "   {"
        "      U               Euler;"
        "      T               EulerExplicit;"
        "   }"
        "}"
        "\n"
        "subModels"
        "{"
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
            "   LiquidPropertiesModel   OFLiquidProperties;"
            "   OFLiquidPropertiesCoeffs"
            "   {"
            "       fluidName   C7H16;"
            "   }"
            "\n"
            "   standardWallInteractionCoeffs"
            "   {"
            "       type            rebound;"
            "   }"
            "\n"
            "   nearestNeighborCoeffs\n"
            "   {\n"
            "       useSPProps                 false;  // false: fLES is used; true: zSP is used\n"
            "       lambda_x                     1.0;"
            "       lambda_f                     0.0;"
            "       lambda_T                     1.0;"
            "       ri                    2.5856e-04;"
            "       fm                        0.0072;"
            "       Tm                          40.0;"
            "   }"
            "   dropletComponent NC7H16;\n"
        "}"
    );

    IOdictionary sprayCloudProperties
    (
        IOobject
        (
            "sprayCloud_nHeptaneProperties",
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        is
    );

    sprayCloudProperties.regIOobject::write();
}


void nHeptaneTest::writeChemistryProperties(Foam::Time& runTime)
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


void nHeptaneTest::writeVolFields
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
        dimensioned<vector>("u",dimVelocity,UDefault)
    );
    U.write();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

