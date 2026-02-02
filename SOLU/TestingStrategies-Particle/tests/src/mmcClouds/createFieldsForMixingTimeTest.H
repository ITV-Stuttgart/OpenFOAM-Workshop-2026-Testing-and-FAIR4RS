
// **********************************************************************
//                  Create fields for mixingModel tests
// **********************************************************************
Info << "Create fields for mixing model tests" << endl;


// Create volVectorField U
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
    dimensioned<vector>("u",dimVelocity,vector(0,0,0))
);
U.write();

// Pressure field
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
    dimensionedScalar("f",dimPressure,1E+5)
);
p.write(); // Needs to be written so basicMMCVar can find this file

// Pressure field
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
    dimensionedScalar("T",dimTemperature,300)
);
T.write();


// Mixture fraction field
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
    dimensionedScalar("f",dimless,0.2)
);
f.write(); // Needs to be written so basicMMCVar can find this file

// Default Y
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
    dimensionedScalar("Ydefault",dimless,0)
);
Ydefault.write();

// Species CH4
volScalarField CH4
(
    IOobject
    (
        "CH4",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("CH4",dimless,0.2)
);
CH4.write();


// Species O2
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
    dimensionedScalar("O2",dimless,0.8)
);
O2.write();

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
    mesh,
    dimensionedScalar("rho0",dimDensity,1)
);
rho.write();



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
dpdt.write();

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
    dimensionedScalar("rho0",dimless,1E-6)
);
D.write();

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
    dimensionedScalar("rho0",dimless,5E-6)
);
Dt.write();

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
DEff.write();

// DeltaE is the eulerian field size
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
DeltaE.ref() = pow(mesh.V(),1.0/3.0);


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

Info << "done"<<endl;
