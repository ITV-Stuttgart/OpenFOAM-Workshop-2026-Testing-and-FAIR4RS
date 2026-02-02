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
    Test the filterDNSField class.

    Uses a function f(x,y,z) that is filtered to check if the filtering function
    works in single core and parallel. 

    Tests a simple function,
    \f[
        \overline{\phi} = \frac{1}{V}\int_{z_0}^{z_1}\int_{y_0}^{y_1}\int_{x_0}^{x_1} x \rd{x}\rd{y}\rd{z}
    \f]
    which is,
    \f[
        \overline{\phi} = \frac{1}{V}\left[0.5*(x_1^2-x_0^2)(y_1-y_0)(z_0-z_1)\right]
    \f]

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2024
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "globalFoamArgs.H"

#include "filterDNSField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace filterFieldsTest
{
    void getCellExtend
    (
        scalar& x0,
        scalar& x1,
        scalar& y0,
        scalar& y1,
        scalar& z0,
        scalar& z1,
        const pointField& points
    );
}


TEST_CASE("filterDNSField Test","[filterFields][filterFields-parallel]")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object
    
    // Construct normal mesh
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

    // Create an own time and object registry for the filter mesh
    Time filterTime(Time::controlDictName,args);    

    // Construct filter mesh
    fvMesh filterMesh
    (
        IOobject
        (
            "filterMesh",
            filterTime.timeName(),
            filterTime,
            IOobject::MUST_READ
        )
    );

    filterDNSField filter(mesh,filterMesh);

    SECTION("Mean Field Test")
    {
        // Create a field on the DNS mesh
        volScalarField f
        (
            IOobject
            (
                "f",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar("f",dimless,0)
        );

        // Populate the field with the x position of the cell center.
        const auto& pos = mesh.C();
        forAll(f,i)
        {
            f[i] = pos[i].x();
        }

        tmp<volScalarField> fFiltered = filter.mean(f);

        // First check that the field has the correct size
        REQUIRE(fFiltered().size() == filterMesh.C().size());

        // Loop over the filter mesh
        // =========================

        // First get the filter mesh cells
        const cellList& cells = filterMesh.cells();

        forAll(filterMesh.C(),i)
        {
            // Compute the extension of the cell
            scalar x0,x1,y0,y1,z0,z1;
            pointField cellPoints = cells[i].points(filterMesh.faces(),filterMesh.points());
            filterFieldsTest::getCellExtend(x0,x1,y0,y1,z0,z1,cellPoints);

            // Integrate the function in space
            scalar integratedValue = 0.5*(std::pow(x1,2)-std::pow(x0,2))*(y1-y0)*(z1-z0)/filterMesh.V()[i];
            REQUIRE_THAT
            (
                integratedValue,
                Catch::Matchers::WithinRel(fFiltered()[i],0.01)
            );

        }

        // Rerun on the same field to check the updateField function
        filter.updateField(f,fFiltered.ref());

        // First check that the field has the correct size
        REQUIRE(fFiltered().size() == filterMesh.C().size());

        // Loop over the filter mesh
        // =========================
        
        INFO("Update field test");
        forAll(filterMesh.C(),i)
        {
            // Compute the extension of the cell
            scalar x0,x1,y0,y1,z0,z1;
            pointField cellPoints = cells[i].points(filterMesh.faces(),filterMesh.points());
            filterFieldsTest::getCellExtend(x0,x1,y0,y1,z0,z1,cellPoints);

            // Integrate the function in space
            scalar integratedValue = 0.5*(std::pow(x1,2)-std::pow(x0,2))*(y1-y0)*(z1-z0)/filterMesh.V()[i];
            REQUIRE_THAT
            (
                integratedValue,
                Catch::Matchers::WithinRel(fFiltered()[i],0.01)
            );
        }
    }

    SECTION("Favre average test")
    {
        // Create a field on the DNS mesh
        volScalarField f
        (
            IOobject
            (
                "f",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar("f",dimless,0)
        );

        // Populate the field with the x position of the cell center.
        const auto& pos = mesh.C();
        forAll(f,i)
        {
            f[i] = pos[i].x();
        }

        volScalarField rho
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar("rho",dimless,0.5)
        );

        tmp<volScalarField> fFiltered = filter.FavreAverage(f,rho);

        // First check that the field has the correct size
        REQUIRE(fFiltered().size() == filterMesh.C().size());

        // Loop over the filter mesh
        // =========================

        // First get the filter mesh cells
        const cellList& cells = filterMesh.cells();

        forAll(filterMesh.C(),i)
        {
            // Compute the extension of the cell
            scalar x0,x1,y0,y1,z0,z1;
            pointField cellPoints = cells[i].points(filterMesh.faces(),filterMesh.points());
            filterFieldsTest::getCellExtend(x0,x1,y0,y1,z0,z1,cellPoints);

            // Integrate the function in space
            scalar integratedValue = 0.5*(std::pow(x1,2)-std::pow(x0,2))*(y1-y0)*(z1-z0)/filterMesh.V()[i];
            REQUIRE_THAT
            (
                integratedValue,
                Catch::Matchers::WithinRel(fFiltered()[i],0.01)
            );

        }
    }


    SECTION("Favre average vector test")
    {
        // Create a field on the DNS mesh
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedVector("U",dimless,vector::zero)
        );

        // Populate the field with the x position of the cell center.
        const auto& pos = mesh.C();
        forAll(U,i)
        {
            U[i] = pos[i];
        }

        volScalarField rho
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            mesh,
            dimensionedScalar("rho",dimless,0.5)
        );

        tmp<volVectorField> UFiltered = filter.FavreAverage(U,rho);

        // First check that the field has the correct size
        REQUIRE(UFiltered().size() == filterMesh.C().size());

        // Loop over the filter mesh
        // =========================

        // First get the filter mesh cells
        const cellList& cells = filterMesh.cells();

        const auto& filterCellPos = filterMesh.C();

        forAll(filterMesh.C(),i)
        {
            // Compute the extension of the cell
            scalar x0,x1,y0,y1,z0,z1;
            pointField cellPoints = cells[i].points(filterMesh.faces(),filterMesh.points());
            filterFieldsTest::getCellExtend(x0,x1,y0,y1,z0,z1,cellPoints);

            // Integrate the function in space
            scalar integratedValueX = 0.5*(std::pow(x1,2)-std::pow(x0,2))*(y1-y0)*(z1-z0)/filterMesh.V()[i];
            scalar integratedValueY = 0.5*(std::pow(y1,2)-std::pow(y0,2))*(x1-x0)*(z1-z0)/filterMesh.V()[i];
            scalar integratedValueZ = 0.5*(std::pow(z1,2)-std::pow(z0,2))*(y1-y0)*(x1-x0)/filterMesh.V()[i];
            INFO("Check cell at position "
                << filterCellPos[i].x() << ", "
                << filterCellPos[i].y() << ", "
                << filterCellPos[i].z());
            REQUIRE_THAT
            (
                integratedValueX,
                Catch::Matchers::WithinRel(UFiltered()[i].x(),0.01)
            );

            REQUIRE_THAT
            (
                integratedValueY,
                Catch::Matchers::WithinRel(UFiltered()[i].y(),0.01)
            );

            REQUIRE_THAT
            (
                integratedValueZ,
                Catch::Matchers::WithinRel(UFiltered()[i].z(),0.01)
            );

        }
    }


    SECTION("mapFilterFieldToDNS Test")
    {
        // Create a filter field on the filter mesh
        volScalarField fFiltered
        (
            IOobject
            (
                "fFiltered",
                filterTime.timeName(),
                filterTime,
                IOobject::NO_READ
            ),
            filterMesh,
            dimensionedScalar("f",dimless,0)
        );

        // Populate the field with the x position of the cell center.
        const auto& pos = filterMesh.C();
        forAll(fFiltered,i)
        {
            fFiltered[i] = pos[i].x();
        }

        tmp<volScalarField> fDNS = filter.mapFilterFieldToDNS(fFiltered);

        // First check that the field has the correct size
        REQUIRE(fDNS().size() == mesh.C().size());

        // Do not write because otherwise the decompose in the runAll does not 
        // work without prior cleaning of the 0/ directory
        // fDNS().write();
        // fFiltered.write();

        // // Loop over the filter mesh
        // // =========================

        // // First get the filter mesh cells
        // const cellList& cells = filterMesh.cells();

        // forAll(filterMesh.C(),i)
        // {
        //     // Compute the extension of the cell
        //     scalar x0,x1,y0,y1,z0,z1;
        //     pointField cellPoints = cells[i].points(filterMesh.faces(),filterMesh.points());
        //     filterFieldsTest::getCellExtend(x0,x1,y0,y1,z0,z1,cellPoints);

        //     // Integrate the function in space
        //     scalar integratedValue = 0.5*(std::pow(x1,2)-std::pow(x0,2))*(y1-y0)*(z1-z0)/filterMesh.V()[i];
        //     REQUIRE_THAT
        //     (
        //         integratedValue,
        //         Catch::Matchers::WithinRel(fFiltered()[i],0.01)
        //     );

        // }
    }

}



void filterFieldsTest::getCellExtend
(
    Foam::scalar& x0,
    Foam::scalar& x1,
    Foam::scalar& y0,
    Foam::scalar& y1,
    Foam::scalar& z0,
    Foam::scalar& z1,
    const Foam::pointField& points
)
{
    x0 =  1E+10;
    x1 = -1E+10;
    y0 =  1E+10;
    y1 = -1E+10;
    z0 =  1E+10;
    z1 = -1E+10;
    for (const auto& p : points)
    {
        if (p.x() > x1)
            x1 = p.x();
        if (p.x() < x0)
            x0 = p.x();

        if (p.y() > y1)
            y1 = p.y();
        if (p.y() < y0)
            y0 = p.y();

        if (p.z() > z1)
            z1 = p.z();
        if (p.z() < z0)
            z0 = p.z();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

