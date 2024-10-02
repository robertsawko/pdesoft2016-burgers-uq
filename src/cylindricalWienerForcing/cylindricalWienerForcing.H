/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Class
    Foam::fv::cylindricalWienerForcing

Group
    grpFvOptionsSources

Description
    Applies a random forcing based on cylindrical Wiener process.

    Sources applied to:
    \verbatim
      U         | Velocity                                   [m/s]
    \endverbatim

    Required fields:
    \verbatim
      U         | Velocity                                   [m/s]
    \endverbatim

Usage
    Minimal example by using \c constant/fvOptions:
    \verbatim
    cylindricalWienerForcing1
    {
        // Mandatory entries (unmodifiable)
        type                 cylindricalWienerForcing;

        // Mandatory entries (unmodifiable)
        Ubar             (1 0 0);
        fields           (<fieldName>);

        // Optional entries (unmodifiable)
        relaxation       1.0;

        // Conditional mandatory entries (unmodifiable)

            // when <timePath>/uniform/<name>Properties file exists
            gradient    <vectorField>; // reading from the aforementioned file

        // Mandatory/Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                       | Type   | Reqd  | Dflt
      type      | Type name: cylindricalWienerForcing      | word   | yes   | -
      Ubar      | Desired mean velocity             | vector | yes   | -
      fields    | Name of operand velocity field    | word   | yes   | -
      relaxation | Relaxation factor                | scalar | no    | 1
      gradient  | Initial pressure gradient field   | vectorField | cndtnl | -
    \endtable

    The inherited entries are elaborated in:
     - \link fvOption.H \endlink
     - \link cellSetOption.H \endlink

Note
  - Currently only handles kinematic pressure (incompressible solvers).

SourceFiles
    cylindricalWienerForcing.C

\*---------------------------------------------------------------------------*/

#ifndef cylindricalWienerForcing_H
#define cylindricalWienerForcing_H

#include "autoPtr.H"
#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
               Class cylindricalWienerForcing Declaration
\*---------------------------------------------------------------------------*/

class cylindricalWienerForcing
:
    public fv::cellSetOption
{
protected:

    // Protected Data

        //- Pressure gradient before correction
        scalar sigma_;

        //- Flow direction
        vector flowDir_;

    // Protected Member Functions

        //- Write the pressure gradient to file (for restarts etc)
        void writeProps(const scalar gradP) const;

        //- Correct driving force for a constant mass flow rate
        void update(fvMatrix<vector>& eqn);


public:

    //- Runtime type information
    TypeName("cylindricalWienerForcing");


    // Constructors

        //- Construct from explicit source name and mesh
        cylindricalWienerForcing
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );

        //- No copy construct
        cylindricalWienerForcing(const cylindricalWienerForcing&) = delete;

        //- No copy assignment
        void operator=(const cylindricalWienerForcing&) = delete;


    //- Destructor
    virtual ~cylindricalWienerForcing() = default;


    // Member Functions

        // Evaluation

            //- Correct the pressure gradient
            virtual void correct(volVectorField& U);

            //- Add explicit contribution to momentum equation
            virtual void addSup
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Add explicit contribution to compressible momentum equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<vector>& eqn,
                const label fieldi
            );

            //- Set 1/A coefficient
            virtual void constrain
            (
                fvMatrix<vector>& eqn,
                const label fieldi
            );


        // IO

            //- Read source dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
