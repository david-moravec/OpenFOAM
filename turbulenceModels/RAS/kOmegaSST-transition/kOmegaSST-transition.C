/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2018 OpenFOAM Foundation
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

#include "kOmegaSST-transition.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegeSST-trasition<BasicTurbulenceModel>::P_gamma
(
    const tmp<volScalarField>& S
    const tmp<volScalarField>& Re_v
    const volScalarField& Rt
)
{
    tmp<volScalarField> p_gamma = F_length() * ca1_ * this->rho * S 
                                * sqrt(gamma_ * F_onset(Re_v, Rt))
                                * (1 - ce1_ * gamma_)

    return p_gamma
}


template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSt-transition<BasicTurbulenceModel>::F_onset
(
    const tmp<volScalarField>& Re_v
    const volScalarField& Rt
)
{
    tmp<volScalarField> F_onset1 = Re_v / (2.193 * this->Re_theta_);
    tmp<volScalarField> F_onset2 = min(max(F_onset1, sqr(sqr(F_onset1))), 2);
    Rt = this->rho_ * this->k_ / (this->nu_ * this->omega_); 
    tmp<volScalarField> F_onset3 = max(1 - pow(Rt/2.5, 3), 0);
    tmp<volScalarField> F_onset = max(F_onset2 - F_onset3, 0);
    return F_onset
}

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSt-transition<BasicTurbulenceModel>::F_onset3

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSt-transition<BasicTurbulenceModel>::F_onset4

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSt-transition<BasicTurbulenceModel>::F_onset2

template<class BasicTurbulenceModel>
kOmegaSST-transition<BasicTurbulenceModel>::kOmegaSST-transition
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    Foam::kOmegaSST-transition
    <
        eddyViscosity<RASModel<BasicTurbulenceModel>>,
        BasicTurbulenceModel
    >
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
