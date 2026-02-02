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
    Foam::EvapModel

Description
    Base class for evaporation models of droplets

Usage
    Minimal example by using \c constant/\<CloudProperties\>:
    \verbatim
    subModels
    {
        EvapModel    none;
    }
    \endverbatim

SourceFiles
    EvapModel.C
    EvapModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef EvapModel_H
#define EvapModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class EvapModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EvapModel
:
    public CloudSubModelBase<CloudType>
{
public:

    // Convienence typedef
    typedef typename CloudType::parcelType::trackingData parcelTrackingData;

    //- Runtime type information
    TypeName("EvapModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        EvapModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Generated Methods

        //- No copy assignment
        void operator=(const EvapModel<CloudType>&) = delete;


    // Constructors

        //- Construct null from owner
        EvapModel(CloudType& owner);

        //- Construct from dictionary
        EvapModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Copy construct
        EvapModel(const EvapModel<CloudType>& satModel);

        //- Construct and return a clone
        virtual autoPtr<EvapModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~EvapModel() = default;


    //- Selector
    static autoPtr<EvapModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

        //- Calculate the mass and heat transfer due to evaporation
        virtual void mDotAndQDot
        (
            const parcelTrackingData& td,   // Tracking data of the parcel
            const scalar& Td,       // Droplet surface temperature
            const vector& Ud,       // Droplet velocity
            const scalar& d,        // Droplet diameter
            const scalar& mD,       // Droplet mass
            const scalar& CpL,      // Droplet heat capacity
            scalar& mDot,
            scalar& QDot
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeEvapModel(CloudType)                                               \
                                                                               \
    typedef Foam::CloudType::EvapDropletCloudType EvapDropletCloudType;        \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::EvapModel<EvapDropletCloudType>,                                 \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            EvapModel<EvapDropletCloudType>,                                   \
            dictionary                                                         \
        );                                                                     \
    }


#define makeEvapModelType(SS, CloudType)                                       \
                                                                               \
    typedef Foam::CloudType::EvapDropletCloudType EvapDropletCloudType;        \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<EvapDropletCloudType>, 0);    \
                                                                               \
    Foam::EvapModel<EvapDropletCloudType>::                                    \
        adddictionaryConstructorToTable<Foam::SS<EvapDropletCloudType>>        \
            add##SS##CloudType##EvapDropletCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EvapModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
