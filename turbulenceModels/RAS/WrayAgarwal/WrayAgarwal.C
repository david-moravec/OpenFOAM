/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "WrayAgarwal.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
	
template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwal<BasicTurbulenceModel>::F1
(
 	const volScalarField& S
)
{
	volScalarField Max =  max 
	(
	 this->y_ * sqrt(this->R_ * S), 
	 1.5 * this->R_
	);

	volScalarField arg1 = (1 + this->y_ * sqrt(this->R_ * S) / this->nu()) /
			      (1 + sqr( Max  / (20 * this->nu())));
	
	return tanh(pow(arg1, 4));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwal<BasicTurbulenceModel>::WrayAgarwal
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
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class BasicTurbulenceModel>
void WrayAgarwal<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volTensorField gradU = fvc::grad(U);

    const volScalarField S = max
	    (
	 	sqrt(2 * symm(gradU) && symm(gradU)),
		dimensionedScalar("1e-16", inv(dimTime), 1e-16)
	    );

//    tmp<volScalarField> W = sqrt(2 * skew(gradU) && skew(gradU));
//    W.clear();


    const volScalarField F1 = this->F1(S);
    const volScalarField nuEff = this->nuEff(F1);
    const volScalarField Fmi = this->Fmi(this->chi());

    volScalarField C1 = F1 * (C1kOm_ - C1kEps_) + C1kEps_;
    volScalarField CD_RS = fvc::grad(R_) & fvc::grad(S);
    
    tmp<fvScalarMatrix> REqn
    (
        fvm::ddt(alpha, rho, R_)
      + fvm::div(alphaRhoPhi, R_)
      - fvm::laplacian(nuEff, R_)
     ==
        alpha * rho * fvm::SuSp(C1 * S, R_)
      + alpha * rho * F1 * fvm::SuSp(C2kOm_ / S * CD_RS, R_)
      - alpha * rho * (1 - F1) * C2kEps_ * fvm::SuSp(R_ * magSqr(fvc::grad(S)) / sqr(S), R_)
    );

    REqn.ref().relax();
    fvOptions.constrain(REqn.ref());
    solve(REqn);
    fvOptions.correct(R_);
    bound(R_, dimensionedScalar(R_.dimensions(), 0));
    R_.correctBoundaryConditions();

    correctNut(Fmi);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam
