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
    Test the MixingPopeParticle

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared

#include "globalFoamArgs.H"

#include "fvCFD.H"
#include "basicMixingPopeParticle.H"
#include "mmcStatusMessage.H"
#include "mmcVarSet.H"
#include "OFstream.H"
#include "IFstream.H"
#include "unitMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("MixingPopeParticle-Test","[popeParticles][popeParticles-parallel]")
{
    unitMesh uMesh(10);
    const fvMesh& mesh = uMesh.mesh();


    SECTION("Construtctor")
    {
        INFO("Create ItoPopeParticle");
        
        label icell = -1;
        label itetface = -1;
        label itetpt = -1;
            
        mesh.findCellFacePt(mesh.C()[0],icell,itetface,itetpt);
        
        // Create a particle
        basicMixingPopeParticle p
        (
            mesh,
            Foam::barycentric(0.123009,0.504829,0.206735,0.165427),
            icell,
            itetface,
            itetpt    
        );
        
        // Set particle properties 
        // note: dxDeterministic is not written also not set here
        p.dx() = 0.0;
        
        
        SECTION("Particle IO - ASCII")
        {
            if (!Pstream::parRun())
            {
            // Write the particle to a file
            OFstream os("mixingPopeParticle.txt",OFstream::streamFormat::ASCII);
            os << p;
            os << flush;
            // Read from file 
            IFstream is("mixingPopeParticle.txt",OFstream::streamFormat::ASCII);
            basicMixingPopeParticle pNew(mesh, is);
    
            REQUIRE(pNew.dx() == p.dx());
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
                        basicMixingPopeParticle pNew(mesh,fromBuffer);
                        REQUIRE(pNew.dx() == p.dx());
                    }
                }
                
            }
        }
    }
    
    
    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

