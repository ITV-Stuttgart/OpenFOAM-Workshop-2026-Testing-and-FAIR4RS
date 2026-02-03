/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "EvapDropletParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
const std::size_t Foam::EvapDropletParcel<ParcelType>::sizeofFields
(
    2*sizeof(scalar)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::EvapDropletParcel<ParcelType>::EvapDropletParcel
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    ParcelType(mesh, is, readFields, newFormat),
    mass0_(0.0)
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is >> mass0_;
            is >> mDot_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, &mass0_);
            readRawScalar(is, &mDot_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&mass0_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


template<class ParcelType>
template<class CloudType>
void Foam::EvapDropletParcel<ParcelType>::readFields(CloudType& c)
{
    const bool valid = c.size();
    ParcelType::readFields(c);

    IOField<scalar> mass0(c.newIOobject("mass0", IOobject::MUST_READ), valid);
    IOField<scalar> mDot(c.newIOobject("mDot", IOobject::MUST_READ), valid);
 
    label i = 0;
    for (EvapDropletParcel<ParcelType>& p : c)
    {
        p.mass0_ = mass0[i];
        p.mDot_ = mDot[i];

 
        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::EvapDropletParcel<ParcelType>::writeFields(const CloudType& c)
{
    ParcelType::writeFields(c);

    const label np = c.size();
    const bool valid = np;
 
    IOField<scalar> mass0(c.newIOobject("mass0", IOobject::NO_READ), np);
    IOField<scalar> mDot(c.newIOobject("mDot", IOobject::NO_READ), np);
 
    label i = 0;
    for (const EvapDropletParcel<ParcelType>& p : c)
    {
        mass0[i] = p.mass0_;
        mDot[i] = p.mDot_;
 
        ++i;
    }
 
    mass0.write(valid);
    mDot.write(valid);
}


template<class ParcelType>
void Foam::EvapDropletParcel<ParcelType>::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    ParcelType::writeProperties(os, filters, delim, namesOnly);

    #undef  writeProp
    #define writeProp(Name, Value)                                            \
        ParcelType::writeProperty(os, Name, Value, namesOnly, delim, filters)

    writeProp("mass0", mass0_);
    writeProp("mDot",mDot_);

    #undef writeProp
}


template<class ParcelType>
template<class CloudType>
void Foam::EvapDropletParcel<ParcelType>::readObjects
(
    CloudType& c,
    const objectRegistry& obr
)
{
    ParcelType::readObjects(c, obr);

    if (!c.size()) return;


    auto& mass0 = cloud::lookupIOField<scalar>("mass0", obr);
    auto& mDot = cloud::lookupIOField<scalar>("mDot", obr);

    label i = 0;
    for (EvapDropletParcel<ParcelType>& p : c)
    {
        p.mass0_ = mass0[i];
        p.mDot_ = mDot[i];

        ++i;
    }
}


template<class ParcelType>
template<class CloudType>
void Foam::EvapDropletParcel<ParcelType>::writeObjects
(
    const CloudType& c,
    objectRegistry& obr
)
{
    ParcelType::writeObjects(c, obr);

    const label np = c.size();

    if (np > 0)
    {
        auto& mass0 = cloud::createIOField<scalar>("mass0", np, obr);
        auto& mDot = cloud::createIOField<scalar>("mDot", np, obr);

        label i = 0;
        for (const EvapDropletParcel<ParcelType>& p : c)
        {
            mass0[i] = p.mass0_;
            mDot[i] = p.mDot_;

            ++i;
        }
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const EvapDropletParcel<ParcelType>& p
)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const ParcelType&>(p)
            << token::SPACE << p.mass0()
            << token::SPACE << p.mDot();
    }
    else
    {
        os  << static_cast<const ParcelType&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.mass0_),
            EvapDropletParcel<ParcelType>::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
