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

Class
    Foam::RASModels::WrayAgarwalTransition

Description
       WrayAgarwalTransitionCoeffs
        {
            sigmakW      0.72;
            sigmakEps    1.0;
            kappa        0.41;
            C1kOm        0.0829;
            C1kEps       0.1127;
            Cw           8.54;
            Cmu          0.09;
            Cm           8.0;
        }
    \endverbatim

SourceFiles
    WrayAgarwalTransition.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwalTransition_H
#define WrayAgarwalTransition_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class WrayAgarwalTransition Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WrayAgarwalTransition
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
protected:

    // Protected data

        // Model coefficients
			dimensionedScalar Flength_;
            dimensionedScalar C1kOm_;
            dimensionedScalar C1kEps_;
            dimensionedScalar sigmakOm_;
            dimensionedScalar sigmakEps_;
            dimensionedScalar kappa_;
            dimensionedScalar C2kOm_;
            dimensionedScalar C2kEps_;
            dimensionedScalar Cmu_;
            dimensionedScalar Comega_;
            dimensionedScalar Cm_;
            dimensionedScalar Ce2_;
            dimensionedScalar Ca2_;
            dimensionedScalar sigmaGamm_;

        // Fields

            volScalarField R_;
            volScalarField gamma_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;
			const volVectorField& n_;


    // Protected Member Functions

        volScalarField chi() const;
		volScalarField Tu_l(const volScalarField& S) const;
        volScalarField Fmi(const volScalarField& chi) const;
		volScalarField Re_thetac(const volScalarField& S) const;
        volScalarField F1 (const volScalarField& S, const volScalarField& W) const;
		volScalarField F_onset(const volScalarField& Re_v,const volScalarField& S) const;
        volScalarField PR_lim(const volScalarField& Re_v,const volScalarField& W) const;
		volScalarField lambdaThetaL() const;
		volScalarField F_PG() const;

        void correctNut(const volScalarField& Fmi);
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("WrayAgarwalTransition");


    // Constructors

        //- Construct from components
        WrayAgarwalTransition
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );

        //- Disallow default bitwise copy construction
        WrayAgarwalTransition(const WrayAgarwalTransition&) = delete;


    //- Destructor
    virtual ~WrayAgarwalTransition()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for R  and gamma
		tmp<volScalarField> nuEff_gamma() const
			{
				 
				return tmp<volScalarField>(
						new volScalarField("nuEff_gamma", this->nu() + (this->nut_ / this->sigmaGamm_))
						);
			}
        tmp<volScalarField> nuEff_R(const volScalarField& F1) const
			{
				volScalarField sigmaR = F1 * (sigmakOm_ - sigmakEps_) + sigmakEps_;
				return tmp<volScalarField>(
						new volScalarField ("nuEff", sigmaR * this->R_ + this->nu())
						);
			}
    	tmp<volScalarField> Rt() const
			{
				return tmp<volScalarField>
					(
						new volScalarField("Rt", this->R_ / this->nu()) 
					);
			}

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        virtual tmp<volScalarField> k() const;

        virtual tmp<volScalarField> epsilon() const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const WrayAgarwalTransition&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
