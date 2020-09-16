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

#include "WrayAgarwalBaseBase.H"
#include "fvOptions.H"
#include "bound.H"
#include "walBaselDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<BasicTurbulenceModel>::chi() const
{
    return this->R_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<BasicTurbulenceModel>::Fmi
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(this->Cw_));
}



template<class BasicTurbulenceModel>
void WrayAgarwalBase<BasicTurbulenceModel>::correctNut
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
void WrayAgarwalBase<BasicTurbulenceModel>::correctNut()
{
	correctNut(Fmi(this->chi()));
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
WrayAgarwalBase<BasicTurbulenceModel>::WrayAgarwalBase
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
   	 	"C1kw",
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
   		0.1127
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool WrayAgarwalBase<BasicTurbulenceModel>::read()
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


        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
volScalarField WrayAgarwalBase<BasicTurbulenceModel>::nuEff
(
 const volScalarField& F1
) const
{
	volScalarField sigmaR = F1 * (sigmakOm_ - sigmakEps_) + sigmakEps_;
	volScalarField nuEff = sigmaR * this->R_ + this->nu();
	return nuEff;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<BasicTurbulenceModel>::k() const
{
    return volScalarField::New
    (
        "k",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -2, 0, 0), 0)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> WrayAgarwalBase<BasicTurbulenceModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Wray-AgarwalBase model. Returning zero field"
        << endl;

    return volScalarField::New
    (
        "epsilon",
        this->mesh_,
        dimensionedScalar(dimensionSet(0, 2, -3, 0, 0), 0)
    );
}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam
