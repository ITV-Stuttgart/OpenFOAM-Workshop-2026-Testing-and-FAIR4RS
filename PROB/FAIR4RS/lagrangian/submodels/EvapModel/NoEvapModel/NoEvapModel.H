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
    Foam::NoEvapModel

Description
    No flashing model

SourceFiles
    NoEvapModel.C

\*---------------------------------------------------------------------------*/

#ifndef NoEvapModel_H
#define NoEvapModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "EvapModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class NoEvapModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoEvapModel
:
    public EvapModel<CloudType>
{
public:

    // Convienence typedef
    typedef typename CloudType::parcelType::trackingData parcelTrackingData;

    //- Runtime type information
    TypeName("none");

    // Generated Methods

        //- No copy assignment
        void operator=(const NoEvapModel<CloudType>&) = delete;


    // Constructors

        //- Construct from dictionary
        NoEvapModel
        (
            const dictionary& dict,
            CloudType& owner
        );

        //- Copy construct
        NoEvapModel(const NoEvapModel<CloudType>& flashingModel);

        //- Construct and return a clone
        virtual autoPtr<EvapModel<CloudType>> clone() const
        {
            return autoPtr<EvapModel<CloudType>>
            (
                new NoEvapModel<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoEvapModel() = default;



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
    #include "NoEvapModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
