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
    Foam::RASModels::WrayAgarwal

Description
       WrayAgarwalCoeffs
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
    WrayAgarwal.C

\*---------------------------------------------------------------------------*/

#ifndef WrayAgarwal_H
#define WrayAgarwal_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class WrayAgarwal Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class WrayAgarwal
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
protected:

    // Protected data

        // Model coefficients
            dimensionedScalar C1kOm_;
            dimensionedScalar C1kEps_;
            dimensionedScalar C2kOm_;
            dimensionedScalar C2kEps_;
            dimensionedScalar sigmakOm_;
            dimensionedScalar sigmakEps_;
            dimensionedScalar kappa_;
            dimensionedScalar Cmu_;
            dimensionedScalar Cw_;

        // Fields

            volScalarField R_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


    // Protected Member Functions

        tmp<volScalarField> chi() const;

        tmp<volScalarField> Fmi(const volScalarField& chi) const;

        tmp<volScalarField> F1
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
    TypeName("WrayAgarwal");


    // Constructors

        //- Construct from components
        WrayAgarwal
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
        WrayAgarwal(const WrayAgarwal&) = delete;


    //- Destructor
    virtual ~WrayAgarwal()
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
        void operator=(const WrayAgarwal&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#ifdef NoRepository
#   include "WrayAgarwal.C"
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
