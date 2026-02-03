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
\*---------------------------------------------------------------------------*/

#include "NoEvapModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoEvapModel<CloudType>::NoEvapModel
(
    const dictionary& dict,
    CloudType& owner
)
:
    EvapModel<CloudType>(owner)
{}


template<class CloudType>
Foam::NoEvapModel<CloudType>::NoEvapModel
(
    const NoEvapModel<CloudType>& evapModel
)
:
    EvapModel<CloudType>(evapModel.owner_)
{}


template<class CloudType>
void Foam::NoEvapModel<CloudType>::mDotAndQDot
(
    const parcelTrackingData& td,   // Tracking data of the parcel
    const scalar& Td,       // Droplet surface temperature
    const vector& Ud,       // Droplet velocity
    const scalar& d,        // Droplet diameter
    const scalar& mD,       // Droplet mass
    const scalar& CpL,      // Droplet heat capacity
    scalar& mDot,
    scalar& QDot
) const
{
    mDot = 0;
    QDot = 0;
    return;
}

// ************************************************************************* //
