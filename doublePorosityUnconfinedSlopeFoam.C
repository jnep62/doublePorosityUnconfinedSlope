/*-----------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
---------------------------------------------------------------------
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

Application
    doublePorosityUnconfinedSlopeFoam

Description
    Uses PISO algorithm; this solver is a clone of pisoFoam.
    Solves for transient heads in double-porosity unconfined Darcy fields. 
	Version for OpenFOAM 8 (fvOptions changed in v9).
    Boundary conditions: only head values, velocity is not used here
\*-----------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    pisoControl piso(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    // convert to non-dimensional equations
    dimensionedScalar ONE_time("ONE", dimensionSet(0,0,1,0,0,0,0), 1);   
    fileName outputFile("hfx0.txt");    // derivative of hf at x=0, for flux
    OFstream os(outputFile);
    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- PISO loop
        while (piso.correct())
        {
            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                volScalarField otherTerms = (kappa * cosTheta)/ 2. * (pow(hm, 2) - pow(hf, 2));
                dhf = fvc::grad(hf);
                hfx = dhf & vector(1,0,0);

                // Head equation for fracture
                fvScalarMatrix hfEqn
                (
                    ONE_time * fvm::ddt(hf)
                  - cosTheta2 * fvm::laplacian(hf, hf) 
                  - sinTheta * hfx
                  - otherTerms 
                );
                solve(hfEqn);

               // Head equation for matrix
                fvScalarMatrix hmEqn
                (
                    ONE_time * fvm::ddt(hm)
                  + (1. / lamda) * otherTerms
                );
                solve(hmEqn);
            }
        }

        runTime.write();
        scalar hfx0=hfx[0];
        os << "time: " << runTime.timeName() << " hfx@x=0: " << hfx0 << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }
 
    Info<< "End\n" << endl;
    return 0;       
}


// *************************************************************** //
