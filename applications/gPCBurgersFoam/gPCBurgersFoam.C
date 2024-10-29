/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    burgersFoam

Description
    Solves intrusive UQ equations for viscous Burgers' equation.

Reference:
    Xiu, Dongbin. Numerical Methods for Stochastic Computations: A Spectral
    Method Approach (Princeton University Press, 2010).
    https://doi.org/10.1515/9781400835348.
    Section 6.5 page 76-77.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createFvOptions.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            forAll(Uhat, k)
            {
                fvVectorMatrix UhatEqn(
                    fvm::ddt(Uhat[k])
                    -
                    fvm::laplacian(nu, Uhat[k])
                );

                forAll(Uhat, i){
                    forAll(Uhat, j){
                        if(j == k) 
                            UhatEqn += e[i][j][k] * fvm::div(phihat[i], Uhat[j]);
                        else
                            UhatEqn += e[i][j][k] * fvc::div(phihat[i], Uhat[j]);
                    }
                }
                solve(UhatEqn);
                phihat[k] = linearInterpolate(Uhat[k]) & mesh.Sf();
            }

        }
        runTime.write();
    }

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
