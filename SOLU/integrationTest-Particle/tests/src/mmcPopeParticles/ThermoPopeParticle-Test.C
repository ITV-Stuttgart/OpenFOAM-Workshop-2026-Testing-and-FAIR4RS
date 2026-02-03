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
    Test the ThermoPopeParticle

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "basicThermoPopeParticle.H"
#include "mmcStatusMessage.H"
#include "mmcVarSet.H"
#include "OFstream.H"
#include "IFstream.H"
#include "unitMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("ThermoPopeParticle-Test","[popeParticles][popeParticles-parallel]")
{
    unitMesh uMesh(10);
    const fvMesh& mesh = uMesh.mesh();

    SECTION("Construtctor")
    {
        INFO("Create ThermoPopeParticle");
        
        label icell = -1;
        label itetface = -1;
        label itetpt = -1;
            
        mesh.findCellFacePt(mesh.C()[0],icell,itetface,itetpt);
        
        // Create a particle
        basicThermoPopeParticle p
        (
            mesh,
            Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
            icell,
            itetface,
            itetpt    
        );
        
        p.T() = 1.0;
        p.Y() = scalarField(10,2.0);
        p.XiC() = scalarField(10,3.0);
        p.hA() = 4.0;
        p.hEqv() = 5.0;
        p.dpYsource() = scalarField(10,7.0);
        p.dp_hAsource() = 8.0;
        
        REQUIRE(p.T() == 1.0);
        REQUIRE(p.Y().size() == 10);
        REQUIRE(p.XiC().size() == 10);
        REQUIRE(p.hA() == 4.0);
        REQUIRE(p.hEqv() == 5.0);
        REQUIRE(p.dpYsource().size() == 10);
        REQUIRE(p.dp_hAsource() == 8.0);
        
        
        SECTION("Particle IO - ASCII")
        {
            if (!Pstream::parRun())
            {
            // Write the particle to a file
            OFstream os("thermoPopeParticle.txt",OFstream::streamFormat::ASCII);
            os << p;
            os << flush;
            // Read from file 
            IFstream is("thermoPopeParticle.txt",OFstream::streamFormat::ASCII);
            basicThermoPopeParticle pNew(mesh, is);
    
            REQUIRE(p.T() == 1.0);
            REQUIRE(p.Y().size() == 10);
            REQUIRE(p.XiC().size() == 10);
            REQUIRE(p.hA() == 4.0);
            REQUIRE(p.hEqv() == 5.0);
            REQUIRE(p.dpYsource().size() == 10);
            REQUIRE(p.dp_hAsource() == 8.0);
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
                        basicThermoPopeParticle pNew(mesh,fromBuffer);
                        REQUIRE(pNew.T() == 1.0);
                        REQUIRE(pNew.Y().size() == 10);
                        REQUIRE(pNew.XiC().size() == 10);
                        REQUIRE(pNew.hA() == 4.0);
                        REQUIRE(pNew.hEqv() == 5.0);
                        REQUIRE(pNew.dpYsource().size() == 10);
                        REQUIRE(pNew.dp_hAsource() == 8.0);
                    }
                }
                
            }
        }
    }
    
    
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

