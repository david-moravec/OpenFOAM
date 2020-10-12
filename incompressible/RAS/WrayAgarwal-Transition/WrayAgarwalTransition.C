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
#include "WrayAgarwalTransition.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{
#define DIM_SC_GAMM(value) dimensionedScalar("value", dimensionSet(0, 2, -1, 0, 0, 0, 0), value)
#define DIM_SC_R(value) dimensionedScalar("value", dimensionSet(0, 2, -3, 0, 0, 0, 0), value)
#define NO_DIM_SC(value) dimensionedScalar("value", dimensionSet(0, 0, 0, 0, 0, 0, 0), value)
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceModel>::chi() const
{
    return this->R_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceModel>::Fmi 
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(this->Cw_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceModel>::F1 
(
 	const volScalarField& S,
 	const volScalarField& W
) const
{
	volScalarField k = this->nut_ * S / sqrt(Cmu_);
	volScalarField w = S / sqrt(Cmu_);
	volScalarField eta = S * max(1.0 , mag(W / S));

	tmp<volScalarField> arg1 = 0.5 * (this->nu() + this->R_) * sqr(eta) / 
                          max(Cmu_ * k * w, DIM_SC_R(1e-16));
	
	return tanh(pow(arg1, 4));
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceModel>::PR_lim
(
 	const volScalarField& W,
	const volScalarField& Re_v
) const
{
	tmp<volScalarField> PR_lim = 1.5 * W * max(this->gamma_ - DIM_SC_GAMM(0.2), DIM_SC_GAMM(0))
						  * (DIM_SC_GAMM(1.0) - this->gamma_)
						  * min(max(this->nut_ / 2420 - DIM_SC_GAMM(1), DIM_SC_GAMM(0)), DIM_SC_GAMM(3))
						  * max(3 * this->nu()- this->nut_, DIM_SC_GAMM(0));
								  
	return PR_lim;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceModel>::F_onset
(
 	const volScalarField& Re_v,
 	const volScalarField& S
) const
{
	tmp<volScalarField> Re_theta = this->Re_theta(S);
	volScalarField R_T = this->nut_ / this->nu();
		
	tmp<volScalarField> F_onset1 = Re_v / (2.2*Re_theta);
	tmp<volScalarField> F_onset2 = min(F_onset1, 2.0);
	tmp<volScalarField> F_onset3 = max(1 - pow(R_T / 3.5, 3), NO_DIM_SC(0));

	return max(F_onset2 - F_onset3, NO_DIM_SC(0));
}
	
template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceModel>::Re_theta
(
 	const volScalarField& S
) const
{
	volScalarField Tu_l = this->Tu_l(S);
	volScalarField F_PG = this->F_PG();
	tmp<volScalarField> arg = 100 + 1e3*Foam::exp(-1.0*Tu_l*F_PG); 
	return arg;
}


template<class BasicTurbulenceMode>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceMode>::p_param
(
) const
{
	volScalarField volume_diff = this->nu() / sqr(this->y_);
	tmp<volScalarField> arg = - 7.57 * 1e-3 * volume_diff * sqr(y_) / this->nu() 
					   + NO_DIM_SC(0.0128);
	return arg;
}

template<class BasicTurbulenceMode>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceMode>::Tu_l
(
 	const volScalarField& S
) const
{
	tmp<volScalarField> arg = 100 * sqrt(2*R_ / 3.0) / (sqrt(S / 0.3) * this->y_);
	return min(arg, NO_DIM_SC(100));
}

template<class BasicTurbulenceMode>
volScalarField WrayAgarwalTransition<BasicTurbulenceMode>::F_PG
(
) const 
{

	double C_PG1 = 14.68;
	double C_PG2 = -7.34;
	double C_PG3 = 0.0;
	double C_PG1_lim = 1.5;
	double C_PG2_lim = 3.0;

	volScalarField lam_theta_L = this->p_param();
	
	volScalarField arg2 = min(lam_theta_L + NO_DIM_SC(0.0681), NO_DIM_SC(0));
	volScalarField F_PG = lam_theta_L;

	forAll(lam_theta_L, i) {
		if (lam_theta_L[i] >= 0) {
			dimensionedScalar arg = 1 + C_PG1 * lam_theta_L[i];
			F_PG[i] = min(arg.value(), C_PG1_lim);
		}
		else {
			dimensionedScalar arg = 1 + C_PG2 * lam_theta_L[i];
			dimensionedScalar arg2 = min(lam_theta_L[i] + 0.0681, 0);
			F_PG[i] = min(arg.value() + C_PG3 * arg2.value(), C_PG2_lim);
		}
	}
	bound(F_PG, dimensionedScalar("0", lam_theta_L.dimensions(), SMALL));
		
	return F_PG;
}

template<class BasicTurbulenceModel>
void WrayAgarwalTransition<BasicTurbulenceModel>::correctNut
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
void WrayAgarwalTransition<BasicTurbulenceModel>::correctNut()
{
	correctNut(Fmi(this->chi()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwalTransition<BasicTurbulenceModel>::WrayAgarwalTransition
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

   Ce2_ 
   (
		dimensioned<scalar>::lookupOrAddToDict
	(
	 	"Ce2",
		this->coeffDict_,
		50.0
	)
   ),

   Ca2_ 
   (
		dimensioned<scalar>::lookupOrAddToDict
	(
	 	"Ca2",
		this->coeffDict_,
		0.06
	)
   ),
	
   sigmaGamm_ 
   (
		dimensioned<scalar>::lookupOrAddToDict
	(
	 	"sigmaGamm",
		this->coeffDict_,
		1.0
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

   gamma_
   (
   IOobject
   	(
   	 	"gamma",
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
bool WrayAgarwalTransition<BasicTurbulenceModel>::read()
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
volScalarField WrayAgarwalTransition<BasicTurbulenceModel>::nuEff_gamma
(
) const
{
	volScalarField nueff = this->nu() + this->nut_ / this->sigmaGamm_;
	return nueff;
}

template<class BasicTurbulenceModel>
volScalarField WrayAgarwalTransition<BasicTurbulenceModel>::nuEff_R
(
 const volScalarField& F1
) const
{
	volScalarField sigmaR = F1 * (sigmakOm_ - sigmakEps_) + sigmakEps_;
	volScalarField nuEff = sigmaR * this->nut_ + this->nu();
	return nuEff;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceModel>::k() const
{
    return volScalarField::New
    (
        "k",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -2, 0, 0), 0)
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalTransition<BasicTurbulenceModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Wray-Agarwal model. Returning zero field"
        << endl;

    return volScalarField::New
    (
        "epsilon",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -3, 0, 0), 0)
    );
}

template<class BasicTurbulenceModel>
void WrayAgarwalTransition<BasicTurbulenceModel>::correct()
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
	const volScalarField Re_v = rho * sqr(this->y_) * S / this->nu();
//functions for gamma
	const volScalarField PR_lim = this->PR_lim(W, Re_v);
    const volScalarField nuEff_gamma = this->nuEff_gamma();
	const volScalarField F_onset = this->F_onset(Re_v, S);

	const volScalarField F_turb = Foam::exp(0.5 * (-this->nut_/ this->nu()));
	double F_length = 100; 

//functions for R
    const volScalarField F1 = this->F1(S, W);
    const volScalarField nuEff_R = this->nuEff_R(F1);
    const volScalarField Fmi = this->Fmi(this->chi());

    const volScalarField C1 = F1 * (C1kOm_ - C1kEps_) + C1kEps_;
    const volScalarField CD_RS = fvc::grad(R_) & fvc::grad(S);
    const volScalarField SS = magSqr(fvc::grad(S)) / sqr(S);

    tmp<fvScalarMatrix> gammaEqn
    (
        fvm::ddt(alpha, rho, gamma_)
      + fvm::div(alphaRhoPhi, gamma_)
      - fvm::laplacian(alpha * rho * nuEff_gamma, gamma_)
     ==
        alpha * rho * F_onset * F_length * fvm::SuSp(S * (DIM_SC_GAMM(1) - gamma_), gamma_)
      - alpha * rho * F_turb * Ca2_ * fvm::SuSp(W * (Ce2_ * gamma_ - DIM_SC_GAMM(1)) , gamma_)
    );

    gammaEqn.ref().relax();
    fvOptions.constrain(gammaEqn.ref());
    solve(gammaEqn);
    fvOptions.correct(gamma_);
   // bound(gamma_, dimensionedScalar(gamma_.dimensions(), 0));
    gamma_.correctBoundaryConditions();
    
    tmp<fvScalarMatrix> REqn
    (
        fvm::ddt(alpha, rho, R_)
      + fvm::div(alphaRhoPhi, R_)
      - fvm::laplacian(alpha * rho * nuEff_R, R_)
     ==
        alpha * rho * gamma_ * fvm::SuSp(C1 * S, R_)
      + alpha * rho * F1 * fvm::SuSp(C2kOm_ / S * CD_RS, R_)
	  + PR_lim
      - alpha * rho * max(gamma_, 0.1) * (1 - F1) * fvm::SuSp(C2kEps_ * R_ * SS, R_)
    );

    REqn.ref().relax();
    fvOptions.constrain(REqn.ref());
    solve(REqn);
    fvOptions.correct(R_);
    bound(R_, dimensionedScalar(R_.dimensions(), 0));
    R_.correctBoundaryConditions();

    correctNut(Fmi);
}
#undef NO_DIM_SC
#undef DIM_SC_R
#undef DIM_SC_GAMM


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

makeTemplatedTurbulenceModel(transportModelIncompressibleTurbulenceModel, RAS, WrayAgarwalTransition)
