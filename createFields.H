Info<< "Reading field hf\n" << endl;

// head in fracture
volScalarField hf
(
    IOobject
    (
        "hf",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//  head in matrix (read for sloping aquifers, to enable variable initial values)
volScalarField hm
(
    IOobject
    (
        "hm",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// derivative vector of hf (output for debugging)
volVectorField dhf
(
   IOobject
    (
        "dhf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::grad(hf)
);

// x-derivative of hf (output for debugging)
volScalarField hfx
(
   IOobject
    (
        "hfx",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    dhf & vector(1,0,0)
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

Info<< "Reading dimensionless parameter kappa [-]\n" << endl;

dimensionedScalar kappa
(
    transportProperties.lookup("kappa")
);

Info<< "Reading dimensionless parameter lamda [-]\n" << endl;

dimensionedScalar lamda
(
    transportProperties.lookup("lamda")
);

// cos and sin are given explicitly; FOAM's angles are in radians
// cos is given twice, with same value but proper dimensions 
Info<< "Reading dimensionless parameter theta cos [-]\n" << endl;
dimensionedScalar cosTheta
(
    transportProperties.lookup("cosTheta")
);
Info<< "Reading dimensionless parameter theta cos2 [-]\n" << endl;
dimensionedScalar cosTheta2
(
    transportProperties.lookup("cosTheta2")
);

Info<< "Reading dimensionless parameter theta sin [-]\n" << endl;

dimensionedScalar sinTheta
(
    transportProperties.lookup("sinTheta")
);

//#include "createFvOptions.H"


