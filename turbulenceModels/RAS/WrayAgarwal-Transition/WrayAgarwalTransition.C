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
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
volScalarField WrayAgarwalTransition<BasicTurbulenceModel>::chi() const
{
    return this->R_/this->nu();
}

template<class BasicTurbulenceModel>
volScalarField WrayAgarwalTransition<BasicTurbulenceModel>::Fmi 
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(this->Comega_));
}

template<class BasicTurbulenceModel>
volScalarField WrayAgarwalTransition<BasicTurbulenceModel>::F1 
(
 	const volScalarField& S,
 	const volScalarField& W
) const
{
	volScalarField K = this->nut_ * S / sqrt(Cmu_);
	volScalarField om = S / sqrt(Cmu_);
	volScalarField eta = S * max(1.0 , mag(W / S));
	volScalarField bounded_CmuKOm = max(Cmu_ * K * om, 
			                            dimensionedScalar("0", dimensionSet(0, 2, -3, 0, 0, 0, 0), SMALL));

	volScalarField arg1 = 0.5 * (this->nu() + this->R_) * sqr(eta) / bounded_CmuKOm;
	
	return tanh(pow4(arg1));
}

template<class BasicTurbulenceModel>
volScalarField WrayAgarwalTransition<BasicTurbulenceModel>::PR_lim
(
 	const volScalarField& W,
	const volScalarField& Re_v
) const
{
	volScalarField F_Onlim = min(max(Re_v / 2420.0 - 1.0, 0.0), 3.0);
	volScalarField PR_lim = 1.5 * W * max(this->gamma_ - 0.2, 0.0) * (1.0 - this->gamma_) * F_Onlim
						  * max(3.0 * this->nu() - this->nut_, 0.*this->nu());
								  
	return PR_lim;
}

template<class BasicTurbulenceModel>
volScalarField WrayAgarwalTransition<BasicTurbulenceModel>::F_onset
(
 	const volScalarField& Re_v,
 	const volScalarField& S
) const
{
	volScalarField F_onset1 = Re_v / (2.2*Re_thetac(S));
	volScalarField F_onset2 = min(F_onset1, 2.0);
	volScalarField F_onset3 = max(1.0 - pow3(Rt() / 3.5), 0.0);

	return max(F_onset2 - F_onset3, 0.0);
}
	
template<class BasicTurbulenceModel>
volScalarField WrayAgarwalTransition<BasicTurbulenceModel>::Re_thetac
(
 	const volScalarField& S
) const
{
	volScalarField Tu_l = this->Tu_l(S);
	volScalarField F_PG = this->F_PG();
	tmp<volScalarField> arg = 100. + 1000.*exp(-1.0*Tu_l*F_PG); 
	return arg;
}


template<class BasicTurbulenceMode>
volScalarField WrayAgarwalTransition<BasicTurbulenceMode>::lambdaThetaL
(
) const
{
	volScalarField dVdy = fvc::grad(this->U_ & n_) & n_;
	volScalarField arg = -7.57e-3 * dVdy * sqr(y_) / this->nu() + 0.0128;

	arg = min(max(arg, -1.0), 1.0);

	return arg;
}

template<class BasicTurbulenceMode>
volScalarField WrayAgarwalTransition<BasicTurbulenceMode>::Tu_l
(
 	const volScalarField& S
) const
{
	volScalarField arg = 100. * sqrt(2.0*R_ / 3.0) / (sqrt(S / 0.3) * this->y_);
	return min(arg, scalar(100.0));
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

	volScalarField lam_theta_L = this->lambdaThetaL();
	volScalarField F_PG = lam_theta_L;

	forAll(lam_theta_L, i) {
		if (lam_theta_L[i] >= 0) {
			double arg = 1.0 + C_PG1 * lam_theta_L[i];
			F_PG[i] = min(arg, C_PG1_lim);
		}
		else {
			double arg = 1.0 + C_PG2 * lam_theta_L[i];
			double arg2 = min(lam_theta_L[i] + 0.0681, 0.0);
			F_PG[i] = min(arg + C_PG3 * arg2, C_PG2_lim);
		}
	}
	return max(F_PG, scalar(0));
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
    Flength_
    (
       dimensioned<scalar>::lookupOrAddToDict
    (
   	 	"Flength",
    	this->coeffDict_,
    	10.0
	)
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
   		0.1338
   	)
   ),

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

   kappa_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"kappa",
   		this->coeffDict_,
   		0.41	
   	)
   ),

   C2kOm_(C1kOm_ / sqr(kappa_) + sigmakOm_),  
   C2kEps_(C1kEps_ / sqr(kappa_) + sigmakEps_),

   Cmu_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"Cmu",
   		this->coeffDict_,
   		0.09
   	)
   ),

   Comega_
   (
    	dimensioned<scalar>::lookupOrAddToDict
   	(
   	 	"Comega",
   		this->coeffDict_,
   		8.54
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

   y_(wallDist::New(this->mesh_).y()),

   n_(wallDist::New(this->mesh_).n())
{
    //this->modifyCoeff(this->C1kEps_, 0.338);
    //this->modifyCoeff(this->C1kOm_, 0.83);
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
		Flength_.readIfPresent(this->coeffDict());
        C1kOm_.readIfPresent(this->coeffDict());
        C1kEps_.readIfPresent(this->coeffDict());
        sigmakOm_.readIfPresent(this->coeffDict());
        sigmakEps_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());
	    C2kOm_ = C1kOm_ / sqr(kappa_) + sigmakOm_;
	    C2kEps_ = C1kEps_ / sqr(kappa_) + sigmakEps_;
        Cmu_.readIfPresent(this->coeffDict());
        Comega_.readIfPresent(this->coeffDict());
        Cm_.readIfPresent(this->coeffDict());
		Ce2_.readIfPresent(this->coeffDict());
		Ca2_.readIfPresent(this->coeffDict());
		sigmaGamm_.readIfPresent(this->coeffDict());


        return true;
    }
    else
    {
        return false;
    }
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
	volScalarField& R_ = this->R_;
	volScalarField& gamma_ = this->gamma_;

    fv::options& fvOptions(fv::options::New(this->mesh_));
    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volTensorField gradU = fvc::grad(U);
    volScalarField S = sqrt(2 * symm(gradU) && symm(gradU));
    volScalarField W = sqrt(2 * skew(gradU) && skew(gradU));
	bound(S, dimensionedScalar("0", S.dimensions(), SMALL));
	bound(W, dimensionedScalar("0", W.dimensions(), SMALL));

	const volScalarField Re_v = sqr(this->y_) * S / this->nu();
	const volScalarField F_turb = exp(-(Rt()/2));
    const volScalarField F1 = this->F1(S, W);
    const volScalarField C1 = F1 * (C1kOm_ - C1kEps_) + C1kOm_;
    
    const volScalarField RS = fvc::grad(R_) & fvc::grad(S);
    const volScalarField SS = magSqr(fvc::grad(S));
    const volScalarField RR = magSqr(fvc::grad(R_));

    tmp<fvScalarMatrix> REqn
    (
        fvm::ddt(alpha, rho, R_)
      + fvm::div(alphaRhoPhi, R_)
      - fvm::laplacian(alpha * rho * nuEff_R(F1), R_)
     ==
        alpha * rho * gamma_ * C1 * R_*S
      + alpha * rho * min(gamma_, 1.) * F1 * C2kOm_ * RS * R_/S
	  + PR_lim(W, Re_v)
      //- alpha * rho * (1 - F1) * C2kEps_ * min(sqr(R_/S) * SS, Cm_ * RR)
      - alpha * rho * min(gamma_, 1.) *  (1 - F1) * C2kEps_ * fvm::Sp(R_/sqr(S) * SS, R_)
    );

    REqn.ref().relax();
    fvOptions.constrain(REqn.ref());
    solve(REqn);
    fvOptions.correct(R_);
    bound(R_, dimensionedScalar(R_.dimensions(), 0));
    R_.correctBoundaryConditions();

    correctNut();

	const volScalarField Pgamma 
	(
        alpha() * rho() * F_onset(Re_v, S) * Flength_ * S * gamma_
	);
	const volScalarField Egamma
	(
      alpha() * rho() * F_turb * Ca2_ * W * gamma_ 
	);

    tmp<fvScalarMatrix> gammaEqn
    (
        fvm::ddt(alpha, rho, gamma_)
      + fvm::div(alphaRhoPhi, gamma_)
      - fvm::laplacian(alpha * rho * nuEff_gamma(), gamma_)
     ==
		Pgamma - fvm::Sp(Pgamma, gamma_)
	  + Egamma - fvm::Sp(Ce2_ * Egamma, gamma_)
    );

    gammaEqn.ref().relax();
    solve(gammaEqn);
    bound(gamma_, dimensionedScalar(gamma_.dimensions(), 0));
	gamma_.correctBoundaryConditions();

    if (debug && this->runTime_.outputTime()) {
        F_onset(Re_v, S)().write();
        F_PG()().write();
	}
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

makeTemplatedTurbulenceModel(transportModelIncompressibleTurbulenceModel, RAS, WrayAgarwalTransition)
