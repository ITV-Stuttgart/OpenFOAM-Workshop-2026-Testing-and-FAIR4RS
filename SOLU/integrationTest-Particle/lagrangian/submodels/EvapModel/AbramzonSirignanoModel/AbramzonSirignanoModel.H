/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 Jan Wilhelm Gärtner
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
    Foam::AbramzonSirignanoModel

Description
    No flashing model

SourceFiles
    AbramzonSirignanoModel.C

\*---------------------------------------------------------------------------*/

#ifndef AbramzonSirignanoModel_H
#define AbramzonSirignanoModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "EvapModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class AbramzonSirignanoModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class AbramzonSirignanoModel
:
    public EvapModel<CloudType>
{

private:

    // Private Member Functions

        //- Calculate the reference conditions
        void calcReferenceValues
        (
            const CloudType& cloud,
            const scalar pInf,          // Infinity pressure
            const scalar Ts,            // Surface temperature
            const scalar TInf,          // Infinity temperature
            const List<scalar>& Ys,     // Surface species composition
            const List<scalar>& YInf,   // Infinity species composition
            const label jFuel,          // Index of the fuel mass fraction in Y
            scalar& TRef,               // Reference temperature
            scalar& rhoRef,             // Density at reference conditions
            scalar& muRef,              // Viscosity at reference conditions
            scalar& kappaRef,           // heat conductivity at reference condition
            scalar& CpRef,              // heat capacity at reference conditions
            scalar& CpFRef              // heat capacity at reference conditions (pure liquid mass fraction)
        ) const;

        //- Transform a species mass fraction vector Y to the mole fraction X
        List<scalar> YtoX(const List<scalar>& Y, const List<scalar>& MW) const;

        //- Transform a species mass fraction vector Y to the mole fraction X
        List<scalar> XtoY(const List<scalar>& X, const List<scalar>& MW) const;

public:

    // Convienence typedef
    typedef typename CloudType::parcelType::trackingData parcelTrackingData;

    //- Runtime type information
    TypeName("AbramzonSirignano");

    // Generated Methods

        //- No copy assignment
        void operator=(const AbramzonSirignanoModel<CloudType>&) = delete;


    // Constructors

        //- Construct from dictionary
        AbramzonSirignanoModel
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Copy construct
        AbramzonSirignanoModel(const AbramzonSirignanoModel<CloudType>& flashingModel);

        //- Construct and return a clone
        virtual autoPtr<EvapModel<CloudType>> clone() const
        {
            return autoPtr<EvapModel<CloudType>>
            (
                new AbramzonSirignanoModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~AbramzonSirignanoModel() = default;



    // Member Functions

        //- Calculate the mass and heat transfer due to evaporation
        void mDotAndQDot
        (
            const parcelTrackingData& td,   // Tracking data of the parcel
            const scalar& Td,       // Droplet surface temperature
            const vector& Ud,       // Droplet velocity
            const scalar& d,        // Droplet diameter
            const scalar& mD,       // Droplet mass
            const scalar& CpL,      // Droplet heat capacity
            scalar& mDot,
            scalar& QDot
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "AbramzonSirignanoModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
