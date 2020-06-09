/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "rotEqTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotEqTotalPressureFvPatchScalarField::
rotEqTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(p, iF),
    omega_(),
    dp0dr_(),
    rRef_()
{}


Foam::rotEqTotalPressureFvPatchScalarField::
rotEqTotalPressureFvPatchScalarField
(
    const rotEqTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    totalPressureFvPatchScalarField(ptf, p, iF, mapper),
    omega_(ptf.omega_, false),
    dp0dr_(ptf.dp0dr_, false),
    rRef_(ptf.rRef_)
{}


Foam::rotEqTotalPressureFvPatchScalarField::
rotEqTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    totalPressureFvPatchScalarField(p, iF, dict),
    omega_(Function1<vector>::New("omega", dict)),
    dp0dr_(Function1<scalar>::New("dp0dr", dict)),
    rRef_(readScalar(dict.lookup("rRef")))
{}


Foam::rotEqTotalPressureFvPatchScalarField::
rotEqTotalPressureFvPatchScalarField
(
    const rotEqTotalPressureFvPatchScalarField& retppsf
)
:
    totalPressureFvPatchScalarField(retppsf),
    omega_(retppsf.omega_, false),
    dp0dr_(retppsf.dp0dr_, false),
    rRef_(retppsf.rRef_)
{}


Foam::rotEqTotalPressureFvPatchScalarField::
rotEqTotalPressureFvPatchScalarField
(
    const rotEqTotalPressureFvPatchScalarField& retppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    totalPressureFvPatchScalarField(retppsf, iF),
    omega_(retppsf.omega_, false),
    dp0dr_(retppsf.dp0dr_, false),
    rRef_(retppsf.rRef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotEqTotalPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    const vector om = omega_->value(t);
    const scalar r0 = rRef_;
    const scalar pGrad = dp0dr_->value(t);

    vector axisHat = om/mag(om);
    tmp<vectorField> rotationVelocity =
        om ^ (patch().Cf() - axisHat*(axisHat & patch().Cf()));

    tmp<scalarField> r = mag(patch().Cf() ^ axisHat);

    const scalarField pCorr
    (
        pGrad*(r - r0)
    );

    const vectorField Up
    (
        patch().lookupPatchField<volVectorField, vector>(UName())
      + rotationVelocity
    );

    totalPressureFvPatchScalarField::updateCoeffs(p0() + pCorr, Up);
}


void Foam::rotEqTotalPressureFvPatchScalarField::write(Ostream& os) const
{
    totalPressureFvPatchScalarField::write(os);
    writeEntry(os, omega_());
    writeEntry(os, dp0dr_());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        rotEqTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
