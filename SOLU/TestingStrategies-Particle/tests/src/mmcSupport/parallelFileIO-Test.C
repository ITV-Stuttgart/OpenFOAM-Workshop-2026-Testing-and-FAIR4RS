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
    Test the baseParticleDataContainer class

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include <catch2/catch_session.hpp> 
#include <catch2/catch_test_macros.hpp> 
#include <catch2/catch_approx.hpp>          // Approx is needed when floats are compared
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "globalFoamArgs.H"

#include "parallelFileIO.H"
#include "fvCFD.H"
#include "IFstream.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


TEST_CASE("parallelFileIO Test","[mmcSupport][mmcSupport-parallel]")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object
    
    // Each processor writes his processor number and 1 behind

	List<scalarList> list(1);
    list[0].resize(2);
    list[0][0] = Pstream::myProcNo();
    list[0][1] = 1;
    fileName testFile("testFile.dat");
    List<word> header(2);
    header[0] = "Proc";
    header[1] = "-";
    
    parallelIO::writeFile(list,header,testFile);
    
    if (Pstream::master())
    {
        // Check the written file by reading the data
        IFstream ifs("testFile.dat");
        List<List<scalar>> readData(Pstream::nProcs());
        
        // discard first line
        std::string dummy;
        ifs.getLine(dummy);

        for (auto& row : readData)
        {
            // Here we have to fix the read format, as the written file does
            // not include a specifier how it is delimited or structured
            row.resize(2);
            ifs >> row[0];
            ifs >> row[1];
        }
        
        forAll(readData,procI)
        {
            REQUIRE(readData[procI][0] == procI);
            REQUIRE(readData[procI][1] == 1);
        }
    }
}


TEST_CASE("parallelFileIO Binary Test","[mmcSupport][mmcSupport-parallel]")
{
    // This test case will write a binary file with the following file format
    // (here as an example for three processors):
    // time processorID
    // 1    0
    // 1    1
    // 1    2
    // 2    0
    // 2    1
    // 2    2

    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object

    // Each processor writes his processor number and 1 behind
    const double deltaT = 0.1;
    const double startTime = 0.0;
    const label nTimeSteps = 10;

	List<scalarList> list(1);
    list[0].resize(2);
    
    fileName testFile("testFileBinary.dat");
    List<word> header(2);
    header[0] = "Time";
    header[1] = "ProcID";

    parallelIO::setStreamPosBasedOnTime(startTime,0,testFile);

    if (Foam::exists("testFileBinary.dat"))
    {
        FatalError 
            << "Delete testFileBinary.dat prior to running the test"
            << exit(FatalError);
    }


    // Write n time steps
    for (label n=0; n < nTimeSteps; n++)
    {
        list[0][0] = startTime+n*deltaT;

        list[0][1] = Pstream::myProcNo();
        parallelIO::writeFileBinary(list,header,testFile);
    }
    
    // Check that it is read correctly
    if (Pstream::master())
    {
        // Check the written file by reading the data
        IFstream ifs("testFileBinary.dat");
        List<List<scalar>> readData;
        
        // discard first line
        std::string dummy;
        ifs.getLine(dummy);
        std::stringstream iss(dummy);

        // Split the string at seperator
        std::string token;
        label nEntries = 0;
        while(std::getline(iss,token,'\t'))
            nEntries++;

        while (ifs.good())
        {
            // Now read binary
            List<scalar> row(nEntries);
            ifs.readRaw(reinterpret_cast<char*>(row.data()),sizeof(scalar)*nEntries);
            if (ifs.good())
                readData.append(row);
        }
        
        REQUIRE(readData.size() == nTimeSteps*Pstream::nProcs());

        for (label n=0; n < nTimeSteps; n++)
        {
            for (label procI=n*Pstream::nProcs(); procI < n*Pstream::nProcs()+Pstream::nProcs(); procI++)
            {
                REQUIRE(readData[procI][0] == startTime+n*deltaT);
                REQUIRE(readData[procI][1] == procI-n*Pstream::nProcs());
            }
        }
    }
}


TEST_CASE("parallelFileIO Binary Test Restart","[mmcSupport]")
{
    // This test needs to be executed after "parallelFileIO Binary Test",
    // and in a seperate execution to reset the static parallelIOFile object!
    //
    //
    // This test case will extend the binary file and add the number 10 to the 
    // processor ID
    // (here as an example for three processors):
    // time processorID
    // 1    10
    // 1    11
    // 1    12
    // 2    10
    // 2    11
    // 2    12

    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object


    //Check if file exists and throw error if it does not exist
    if (!Foam::exists("testFileBinary.dat"))
    {
        FatalError 
            << "File testFileBinary.dat must be generated by the test case "
            << "\"parallelFileIO Binary Test\" prior to this run"
            << exit(FatalError);
    }

    // Each processor writes his processor number and 1 behind
    const double deltaT = 0.1;
    const double startTime = 0.5;
    const label nTimeSteps = 10;

	List<scalarList> list(1);
    list[0].resize(2);
    
    fileName testFile("testFileBinary.dat");
    List<word> header(2);
    header[0] = "Time";
    header[1] = "ProcID";

    // Set the stream position based on the time 
    parallelIO::setStreamPosBasedOnTime(startTime,0,testFile);


    // Write n time steps
    for (label n=0; n < nTimeSteps; n++)
    {
        list[0][0] = startTime+n*deltaT;

        list[0][1] = 10+Pstream::myProcNo();
        parallelIO::writeFileBinary(list,header,testFile);
    }
    
    // Check that it is read correctly
    if (Pstream::master())
    {
        // Check the written file by reading the data
        IFstream ifs("testFileBinary.dat");
        List<List<scalar>> readData;
        
        // discard first line
        std::string dummy;
        ifs.getLine(dummy);
        std::stringstream iss(dummy);

        // Split the string at seperator
        std::string token;
        label nEntries = 0;
        while(std::getline(iss,token,'\t'))
            nEntries++;

        while (ifs.good())
        {
            // Now read binary
            List<scalar> row(nEntries);
            ifs.readRaw(reinterpret_cast<char*>(row.data()),sizeof(scalar)*nEntries);
            if (ifs.good())
                readData.append(row);
        }

        // Check the old, unmodified content
        for (label n=0; n < startTime/deltaT; n++)
        {
            for (label rowI=n*Pstream::nProcs(); rowI < n*Pstream::nProcs()+Pstream::nProcs(); rowI++)
            {
                REQUIRE(readData[rowI][0] == n*deltaT);
                REQUIRE(readData[rowI][1] == rowI-n*Pstream::nProcs());
            }
        }

        for (label n=startTime/deltaT; n < nTimeSteps; n++)
        {
            for (label rowI=n*Pstream::nProcs(); rowI < n*Pstream::nProcs()+Pstream::nProcs(); rowI++)
            {
                REQUIRE_THAT(readData[rowI][0], Catch::Matchers::WithinULP(n*deltaT,1));
                REQUIRE(readData[rowI][1] == 10+rowI-n*Pstream::nProcs());
            }
        }
    }
}









// Following test case can be used to read the data and check 
// manually if it is correct
TEST_CASE("parallelFileIO Manual Check")
{
    // Replace setRootCase.H for Catch2   
    Foam::argList& args = getFoamArgs();
    #include "createTime.H"        // create the time object
        
    if (Pstream::master())
    {
        // Check the written file by reading the data
        IFstream ifs("testFileBinary.dat");
        List<List<scalar>> readData;
        
        // discard first line
        std::string dummy;
        ifs.getLine(dummy);
        std::stringstream iss(dummy);

        // Split the string at seperator
        std::string token;
        label nEntries = 0;
        while(std::getline(iss,token,'\t'))
            nEntries++;

        while(ifs.good())
        {
            // Now read binary
            List<scalar> row(nEntries);
            ifs.readRaw(reinterpret_cast<char*>(row.data()),sizeof(scalar)*nEntries);
            if (ifs.good())
            {
                Info << "row: " << row << endl;
                readData.append(row);
            }
        }
        Info << "nRows: "<<readData.size()<<endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

