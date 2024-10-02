/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "cylindricalWienerForcing.H"
#include "fvMatrices.H"
#include "DimensionedField.H"
#include "IFstream.H"
#include "addToRunTimeSelectionTable.H"
#include "volFieldsFwd.H"
#include "FieldOps.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(cylindricalWienerForcing, 0);
    addToRunTimeSelectionTable(option, cylindricalWienerForcing, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::cylindricalWienerForcing::writeProps
(
    const scalar gradP
) const
{
    // Only write on output time
    if (mesh_.time().writeTime())
    {
        IOdictionary propsDict
        (
            IOobject
            (
                name_ + "Properties",
                mesh_.time().timeName(),
                "uniform",
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            )
        );
        propsDict.add("gradient", gradP);
        propsDict.regIOobject::write();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::cylindricalWienerForcing::cylindricalWienerForcing
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    sigma_(coeffs_.getOrDefault<scalar>("sigma", 1.0)),
    flowDir_(coeffs_.getOrDefault<vector>("flowDir", vector(1.0, 0.0, 0.0)))
{
    /*coeffs_.readEntry("fields", fieldNames_);

    if (fieldNames_.size() != 1)
    {
        FatalErrorInFunction
            << "settings are:" << fieldNames_ << exit(FatalError);
    }
    */

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::cylindricalWienerForcing::correct(volVectorField& U)
{
}


void Foam::fv::cylindricalWienerForcing::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    volVectorField::Internal Su
    (
        IOobject
        (
            name_ + fieldNames_[fieldi] + "Sup",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector(eqn.dimensions()/dimVolume, Zero)
    );

    volScalarField dW (
        IOobject
        (
         "dW",
         mesh_.time().timeName(),
         mesh_,
         IOobject::NO_READ,
         IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimAcceleration, Zero)
    );
    Random::gaussianGeneratorOp<scalar> gen(mesh_.time().timeIndex());
    FieldOps::assign(dW, dW, gen);
    dW *= sigma_/0.01;


    UIndirectList<vector>(Su, cells_) = flowDir_*dW;

    eqn += Su;
}


void Foam::fv::cylindricalWienerForcing::addSup
(
    const volScalarField& rho,
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    this->addSup(eqn, fieldi);
}


void Foam::fv::cylindricalWienerForcing::constrain
(
    fvMatrix<vector>& eqn,
    const label
)
{
}


bool Foam::fv::cylindricalWienerForcing::read(const dictionary& dict)
{
    NotImplemented;

    return false;
}


// ************************************************************************* //
