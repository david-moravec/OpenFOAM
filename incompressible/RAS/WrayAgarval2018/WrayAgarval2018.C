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

#include "WrayAgarval2018.H"
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
tmp<volScalarField> WrayAgarval2018<BasicTurbulenceModel>::chi() const
{
    return this->R_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarval2018<BasicTurbulenceModel>::Fmi
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(this->Cw_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarval2018<BasicTurbulenceModel>::F1
(
 	const volScalarField& S,
 	const volScalarField& W
)
{
	volScalarField k = this->nut_ * S / sqrt(Cmu_);
	volScalarField w = S / sqrt(Cmu_);
	volScalarField eta = S * max(1.0 , mag(W / S));

	volScalarField arg1 = 0.5 * (this->nu() + this->R_) * sqr(eta) / 
                          max(Cmu_ * k * w,
                              dimensionedScalar("SMALL", 
                                            dimensionSet(0, 2, -3, 0 ,0,0,0),       
                                            SMALL));
	
	return tanh(pow(arg1, 4));
}


/*
template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarval2018<BasicTurbulenceModel>::F1
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

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarval2018<BasicTurbulenceModel>::F1
(
 	const volScalarField& S,
 	const volScalarField& W
)
{
    volScalarField eta = S * max(1.0, mag(W/S));
    volScalarField Om = S / sqrt(Cmu_);
    volScalarField k = this->nut_ * Om;
    
	volScalarField arg1 = (this->nut_ + this->R_) / 2. * sqr(eta) / (Cmu_  * k *  Om); 

	return tanh(pow(arg1, 4));
}

*/
template<class BasicTurbulenceModel>
void WrayAgarval2018<BasicTurbulenceModel>::correctNut
(
 	const volScalarField& Fmi
)
{
	this->nut_ = Fmi * this->R_;
	this->nut_.correctBoundaryConditions();
	fv::options::New(this->mesh_).correct(this->nut_);

	BasicTurbulenceModel::correctNut();
}
	
template<class BasicTurbulenceModel>
void WrayAgarval2018<BasicTurbulenceModel>::correctNut()
{
	correctNut(Fmi(this->chi()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarval2018<BasicTurbulenceModel>::WrayAgarval2018
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

    C1kOm_
    (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"C1kOm",
   		this->coeffDict_,
   		0.0833
   	)
   ),

   C1kEps_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"C1kEps",
   		this->coeffDict_,
   		0.1284
   	)
   ),

   Cm_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"Cm",
   		this->coeffDict_,
   		8.0
   	)
   ),

   C2kOm_(C1kOm_ / sqr(kappa_) + sigmakOm_),  
   C2kEps_(C1kEps_ / sqr(kappa_) + sigmakEps_),

   sigmakOm_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"sigmakOm",
   		this->coeffDict_,
   		0.72	
   	)
   ),

   sigmakEps_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"sigmakEps",
   		this->coeffDict_,
   		1.0	
   	)
   ),

   Cmu_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"Cmu",
   		this->coeffDict_,
   		0.09
   	)
   ),

   Cw_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"Cw",
   		this->coeffDict_,
   		8.54
		)
   ),

   kappa_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"kappa",
   		this->coeffDict_,
   		0.41	
   	)
   ),
	
   R_
   (
   IOobject
   	(
   	 	"R",
   		this->runTime_.timeName(),
   		this->mesh_,
   		IOobject::MUST_READ,
   		IOobject::AUTO_WRITE
   	),
   	this->mesh_
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
bool WrayAgarval2018<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        sigmakOm_.readIfPresent(this->coeffDict());
        sigmakEps_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
        C1kOm_.readIfPresent(this->coeffDict());
        C1kEps_.readIfPresent(this->coeffDict());
	    C2kOm_ = C1kOm_ / sqr(kappa_) + sigmakOm_;
	    C2kEps_ = C1kEps_ / sqr(kappa_) + sigmakEps_;
        Cw_.readIfPresent(this->coeffDict());
        Cmu_.readIfPresent(this->coeffDict());
        Cm_.readIfPresent(this->coeffDict());


        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
volScalarField WrayAgarval2018<BasicTurbulenceModel>::nuEff
(
 const volScalarField& F1
) const
{
	volScalarField sigmaR = F1 * (sigmakOm_ - sigmakEps_) + sigmakEps_;
	volScalarField nuEff = sigmaR * this->R_ + this->nu();
	return nuEff;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarval2018<BasicTurbulenceModel>::k() const
{
    return volScalarField::New
    (
        "k",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -2, 0, 0), 0)
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarval2018<BasicTurbulenceModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Wray-Agarval model. Returning zero field"
        << endl;

    return volScalarField::New
    (
        "epsilon",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -3, 0, 0), 0)
    );
}

template<class BasicTurbulenceModel>
void WrayAgarval2018<BasicTurbulenceModel>::correct()
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

    const volScalarField W = sqrt(2 * skew(gradU) && skew(gradU));


    const volScalarField F1 = this->F1(S, W);
    const volScalarField nuEff = this->nuEff(F1);
    const volScalarField Fmi = this->Fmi(this->chi());

    volScalarField C1 = F1 * (C1kOm_ - C1kEps_) + C1kEps_;
    volScalarField CD_RS = fvc::grad(R_) & fvc::grad(S);
    volScalarField SS_RR = min
                           (
                            C2kEps_ * R_  * R_ * magSqr(fvc::grad(S)) / sqr(S),
                            Cm_ * magSqr(fvc::grad(R_)) 
                            );

    
    tmp<fvScalarMatrix> REqn
    (
        fvm::ddt(alpha, rho, R_)
      + fvm::div(alphaRhoPhi, R_)
      - fvm::laplacian(nuEff, R_)
     ==
        alpha * rho * fvm::SuSp(C1 * S, R_)
      + alpha * rho * F1 * fvm::SuSp(C2kOm_ / S * CD_RS, R_)
      - alpha * rho * (1 - F1) * SS_RR
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

#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"
#include "RASModel.H"
#include "transportModel.H"
#include "incompressibleTurbulenceModel.H"
#include "IncompressibleTurbulenceModel.H"

namespace Foam
{
	typedef IncompressibleTurbulenceModel<transportModel>
		transportModelIncompressibleTurbulenceModel;
	typedef RASModel<transportModelIncompressibleTurbulenceModel>
		RAStransportModelIncompressibleTurbulenceModel;
}

makeTemplatedTurbulenceModel(transportModelIncompressibleTurbulenceModel, RAS, WrayAgarval2018)


// ************************************************************************* //
