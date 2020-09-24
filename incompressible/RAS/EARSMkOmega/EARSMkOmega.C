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

#include "EARSMkOmega.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// ************************ ///new Function// ************************ // 

template<class BasicTurbulenceModel>
tmp<volScalarField> 
EARSMkOmega<BasicTurbulenceModel>::F_TNT
(
 const volScalarField& CDkOmega
) const
{ 
	tmp<volScalarField> argTNT = max
	(
	 	CDkOmega,
		dimensionedScalar(dimensionSet(0,0,-3,0,0,0,0), 0)
	);
	return (this->alphaD_ / omega_ * argTNT);
}

template<class BasicTurbulenceModel>
volScalarField
EARSMkOmega<BasicTurbulenceModel>::N
(
 const volScalarField P1,
 const volScalarField P2,
 const volScalarField A3p
) const {
	volScalarField Np = A3p / 3.;

	forAll(Np, i)
	{
		if (P2[i] >= 0) {
			scalar a = P1[i] + sqrt(P2[i]);
			Np += pow(a, 1./3.) + sign(a) * pow(fabs(a), 1./3.);
		}
		
		else {
			scalar b = pow(P1[i], 2.) - P2[i];
			Np += 2 * pow(b, 1./6.) * cos( 1./3. * acos( P1[i] / sqrt(b) ) );
		}
	};
	return Np;
}



//nonlinear stress
template<class BasicTurbulenceModel>
volTensorField
EARSMkOmega<BasicTurbulenceModel>::correctNonlinearStress(const volTensorField& gradU)
{
	volScalarField tau(1. / (this->Cmu_ * omega_));

	volSymmTensorField S(tau * symm(gradU));
	volTensorField W(tau * skew(gradU));

	volScalarField IIW( tr(W & W) );
	volScalarField IIS( tr(S & S) );
	volScalarField IV( tr(S & W & W) );	

	scalar Cdiff = 2.2;
	scalar Neq = 81./20.;

	volScalarField beta1eq = -6./5. * ( Neq / (pow(Neq, 2.) - 2*IIW) );
	volScalarField A3p =  9./5. + 9./4. * Cdiff * max( 1. + beta1eq * IIS, 0.0 );
	volScalarField P1 = A3p * (pow(A3p, 2)/27. + 9./20. * IIS - 2./3. * IIW);
	volScalarField P2 = sqr(P1) * pow( (sqr(A3p)/9. + 9./10. * IIS + 2./3. * IIW), 3);

	volScalarField N = this->N(P1, P2, A3p);

	volScalarField Q = 5./6. * (sqr(N) - 2 * IIW) * (2 * sqr(N) - IIW);
	volScalarField beta1 = -N * (2* sqr(N) - 7 * IIW) / Q;
	volScalarField beta3 = -12 * IV / (N * Q);
	volScalarField beta4 = -2 * (sqr(N) - 2 * IIW) / Q;
	volScalarField beta6 = -6 * N / Q;
	volScalarField beta9 = 6 / Q;
	
	volScalarField Cmut = -0.5 * (beta1 + IIW * beta6);
    	this->nut_ = Cmut / this->beta_ * this->rho_ * k_/omega_;
    	this->nut_.correctBoundaryConditions();
    	fv::options::New(this->mesh_).correct(this->nut_);

	volTensorField nonlinearStress = beta3 * (W & W - 1./3. * IIW * I) + 
			 		 beta4 * (S & W - W & S) +
					 beta6 * (S & W & W + W & W & S - IIW * S - 2./3. * IV * I) +
					 beta9 * (W & S & W & W - W & W & S & W);
			   
    	BasicTurbulenceModel::correctNut();

	return nonlinearStress;
}

// ****************************************************************** // 
/*
template<class BasicTurbulenceModel>
void EARSMkOmega<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu / this->beta_ * this->rho_ * k_/omega_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}
*/

template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> EARSMkOmega<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> EARSMkOmega<BasicTurbulenceModel>::omegaSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            omega_,
            dimVolume*this->rho_.dimensions()*omega_.dimensions()/dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
EARSMkOmega<BasicTurbulenceModel>::EARSMkOmega
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
    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel>>
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

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            this->coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta",
            this->coeffDict_,
            0.075
        )
    ),
    gamma_    // changed from 0.53 
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma",
            this->coeffDict_,
            0.553
        )
    ),
    alphaK_  // changed
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            this->coeffDict_,
            1.01
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            this->coeffDict_,
            0.5
        )
    ),
    alphaD_ //added because of the new function F_TNT
    (
     	dimensioned<scalar>::lookupOrAddToDict
	(
	    "alphaD",
	    this->coeffDict_,
	    0.52
	)
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    omega_
    (
        IOobject
        (
            IOobject::groupName("omega", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(omega_, this->omegaMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool EARSMkOmega<BasicTurbulenceModel>::read()
{
    if (nonlinearEddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());
        beta_.readIfPresent(this->coeffDict());
        gamma_.readIfPresent(this->coeffDict());
        alphaK_.readIfPresent(this->coeffDict());
        alphaOmega_.readIfPresent(this->coeffDict());
	alphaD_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void EARSMkOmega<BasicTurbulenceModel>::correct()
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
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    nonlinearEddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);

    volTensorField nonlinearStress = this->correctNonlinearStress(tgradU);

    volScalarField::Internal G
    (
        this->GName(),
        nut.v()*(dev(twoSymm(tgradU().v())) && tgradU().v()) - nonlinearStress
    );
    tgradU.clear();

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

// ************************ ///new scalar fields// ************************ // 

    //implemantation of cross difusion term
    volScalarField CDkOmega
    (
	fvc::grad(k_) & fvc::grad(omega_) 
    );
/*
    volScalarField CDkOmega = max
    (
	(this->alphaD_) / omega_  * (fvc::grad(k_) & fvc::grad(omega_)),
        dimensionedScalar(dimensionSet(0,0,-2,0,0,0,0), 0)
    );
*/

    //initializing of F_TNT
    volScalarField F_TNT(this->F_TNT(CDkOmega));

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(alpha, rho, omega_)
      + fvm::div(alphaRhoPhi, omega_)
      - fvm::laplacian(alpha*rho*DomegaEff(), omega_)
     ==
        gamma_*alpha()*rho()*G*omega_()/k_()
      - fvm::SuSp(((2.0/3.0)*gamma_)*alpha()*rho()*divU, omega_)
      - fvm::Sp(beta_*alpha()*rho()*omega_(), omega_)
      + omegaSource()
      + fvOptions(alpha, rho, omega_)
//      +	 alpha * rho * CDkOmega //TNT 
      +	 alpha * rho * F_TNT //TNT 
    );

    omegaEqn.ref().relax();
    fvOptions.constrain(omegaEqn.ref());
    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    solve(omegaEqn);
    fvOptions.correct(omega_);
    bound(omega_, this->omegaMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(Cmu_*alpha()*rho()*omega_(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
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

makeTemplatedTurbulenceModel(transportModelIncompressibleTurbulenceModel, RAS, EARSMkOmega)
// ************************************************************************* //
