/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

template<class ParcelType>
template<class TrackCloudType>
inline Foam::EvapDropletParcel<ParcelType>::trackingData::trackingData
(
    const TrackCloudType& cloud,
    trackPart part
)
:
    ParcelType::trackingData(cloud, part),
    Y_(cloud.thermo().carrier().Y()),
    pInterp_
    (
        interpolation<scalar>::New
        (
            cloud.solution().interpolationSchemes(),
            cloud.p()
        )
    ),
    pc_(Zero),
    Yc_(Y_.size(),0.0)
{
    YiInterp_.resize(Y_.size());
    forAll(Y_,i)
    {   
        YiInterp_[i].reset
        (
            interpolation<scalar>::New
            (
                cloud.solution().interpolationSchemes().template get<word>("Yi"),
                Y_[i]
            )
        );
    }
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::EvapDropletParcel<ParcelType>::trackingData::pInterp() const
{
    return *pInterp_;
}


template<class ParcelType>
inline const Foam::interpolation<Foam::scalar>&
Foam::EvapDropletParcel<ParcelType>::trackingData::YInterp(const label i) const
{
    return *(YiInterp_[i]);
}


template<class ParcelType>
inline Foam::scalar Foam::EvapDropletParcel<ParcelType>::trackingData::pc() const
{
    return pc_;
}


template<class ParcelType>
inline Foam::scalar& Foam::EvapDropletParcel<ParcelType>::trackingData::pc()
{
    return pc_;
}


template<class ParcelType>
inline const Foam::List<Foam::scalar>& Foam::EvapDropletParcel<ParcelType>::trackingData::Yc() const
{
    return Yc_;
}


template<class ParcelType>
inline Foam::List<Foam::scalar>& Foam::EvapDropletParcel<ParcelType>::trackingData::Yc()
{
    return Yc_;
}

// ************************************************************************* //
