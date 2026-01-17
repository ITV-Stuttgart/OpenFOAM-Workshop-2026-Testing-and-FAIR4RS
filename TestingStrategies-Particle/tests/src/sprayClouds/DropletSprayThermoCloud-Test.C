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
    Test the SprayMMCCouplingCloud class

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2024
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "Cloud.H"
#include "basicDropletSprayThermoCloud.H"
#include "basicReactingPopeCloud.H"
#include "basicThermoPopeParticle.H"
#include "mmcVarSet.H"
#include "psiReactionThermo.H"
#include "mmcStatusMessage.H"
#include "dimensionedVector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("DropletSprayThermoCloud-Test","[sprayClouds]")
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
    #include "createFieldsForMMC.H"// create all necessary input fields 

    INFO("Create the MMC cloud");
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


    SECTION("Coded Saturation Pressure")
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
            "          U               semiImplicit 1;"
            "      }"
            "      resetOnStartup   false;"
            "   }"
            "\n"
            "   interpolationSchemes"
            "   {"
            "      rho             cell;"
            "      U               cellPoint;"
            "      muc             cell;"
            "      p               cell;"
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
            "       binaryDiffusionCoefficient  5.765e-06;\n"
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
            "               std::vector<scalar> a {2052.57331394213, -0.121990926653498, -0.00714146172046278, 5.20523562482363e-05, 0.0};  // C2H5OH\n"
            "               return ( a[0]+a[1]*x+a[2]*x*x+a[3]*x*x*x+a[4]*x*x*x*x );\n"
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


        // Create the coupling cloud
        basicDropletSprayThermoCloud sprayCloud
        (
            "sprayCloud",
            rho,
            U,
            g,
            mmcCloud.slgThermo(),
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
        REQUIRE_THAT(sprayCloud.liquidProperties().Df(101325,298,28),Catch::Matchers::WithinRel(1.35E-6,0.1));

        REQUIRE_THAT(sprayCloud.liquidProperties().Lv(101325,320),Catch::Matchers::WithinRel(909809.67,0.1));
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

