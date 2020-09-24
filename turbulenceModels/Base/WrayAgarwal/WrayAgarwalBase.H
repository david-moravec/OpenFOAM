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
    Foam::RASModels::WrayAgarwalBase

Description
       WrayAgarwalBaseCoeffs
        {
            sigmakW      0.72;
            sigmakEps    1.0;
            kappa        0.41;
            C1kOm        0.0829;
            C1kEps       0.1127;
            Cw           8.54;
            Cmu          0.09;
        }
    \endverbatim

SourceFiles
    WrayAgarwalBase.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwalBase_H
#define WrayAgarwalBase_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class WrayAgarwalBase Declaration
\*---------------------------------------------------------------------------*/

template<class TurbulenceModel, class BasicTurbulenceModel>
class WrayAgarwalBase
:
    public TurbulenceModel
{
protected:

    // Protected data

        // Model coefficients
            dimensionedScalar sigmakOm_;
            dimensionedScalar sigmakEps_;
            dimensionedScalar kappa_;
            dimensionedScalar C1kOm_;
            dimensionedScalar C1kEps_;
            dimensionedScalar C2kOm_;
            dimensionedScalar C2kEps_;
            dimensionedScalar Cw_;
            dimensionedScalar Cmu_;

        // Fields

            volScalarField R_;

            //- Wall distance
            //  Note: different to walBasel distance in parent RASModel
            //  which is for near-walBasel cells only


    // Protected Member Functions

        virtual tmp<volScalarField> chi() const;

        virtual tmp<volScalarField> Fmi(const volScalarField& chi) const;

        virtual tmp<volScalarField> F1
        (
            const volScalarField& S
        );

        void correctNut(const volScalarField& fv1);
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("WrayAgarwalBase");


    // Constructors

        //- Construct from components
        WrayAgarwalBase
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
        WrayAgarwalBase(const WrayAgarwalBase&) = delete;


    //- Destructor
    virtual ~WrayAgarwalBase()
    {}


    // Member Functions

        //- Read RASProperties dictionary
        virtual bool read();

        //- Return the effective diffusivity for R 
        volScalarField nuEff
		(
		    const volScalarField& F1
		 ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

        virtual tmp<volScalarField> k() const;

        virtual tmp<volScalarField> epsilon() const;

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const WrayAgarwalBase&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "WrayAgarwalBase.C"
#endif
#endif

// ************************************************************************* //
