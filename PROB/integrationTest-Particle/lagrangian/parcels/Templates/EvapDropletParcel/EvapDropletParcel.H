/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    Foam::EvapDropletParcel

Group
    grpLagrangianIntermediateParcels

Description
    Evaporating liquid parcel class with different evaporation sub models.

SourceFiles
    EvapDropletParcelI.H
    EvapDropletParcel.C
    EvapDropletParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef EvapDropletParcel_H
#define EvapDropletParcel_H

#include "particle.H"
#include "SLGThermo.H"
#include "demandDrivenEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class EvapDropletParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const EvapDropletParcel<ParcelType>&
);


/*---------------------------------------------------------------------------*\
                        Class EvapDropletParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class EvapDropletParcel
:
    public ParcelType
{
    public:

    //- Size in bytes of the fields
    static const std::size_t sizeofFields;
 

    class trackingData
    :
        public ParcelType::trackingData
    {
    private:

        // Private data

            //- Reference of the species mass fraction fields
            const PtrList<volScalarField>& Y_;

            // Interpolators for continuous phase fields

                //- Interpolator for continuous phase pressure field
                autoPtr<interpolation<scalar>> pInterp_;

                //- Reference variable list of field interpolators
                List<autoPtr<interpolation<scalar>>> YiInterp_;

            // Cached continuous phase properties

                //- Pressure [Pa]
                scalar pc_;

                //- Species mass fraction
                List<scalar> Yc_;


    public:

        typedef typename ParcelType::trackingData::trackPart trackPart;

        // Constructors

            //- Construct from components
            template<class TrackCloudType>
            inline trackingData
            (
                const TrackCloudType& cloud,
                trackPart part = ParcelType::trackingData::tpLinearTrack
            );


        // Member functions

            //- Return const access to the interpolator for continuous phase
            //  pressure field
            inline const interpolation<scalar>& pInterp() const;

            //- Return const access to the interpolator for continuous phase
            //  species mass fraction field
            inline const interpolation<scalar>& YInterp(const label i) const;

            //- Return the continuous phase pressure
            inline scalar pc() const;

            //- Access the continuous phase pressure
            inline scalar& pc();

            //- Return the continuous phase species mass fraction
            inline const List<scalar>& Yc() const;

            //- Return the continuous phase species mass fraction
            inline List<scalar>& Yc();
    };


protected:

    // Protected data

        // Parcel properties

            //- Initial mass [kg]
            scalar mass0_;

            //- Store the evaporation rate for post-processing
            scalar mDot_;

public:

    // Static data members

        //- Runtime type information
        TypeName("EvapDropletParcel");

        //- String representation of properties
        AddToPropertyList
        (
            ParcelType,
            " mass0"
          + " mDot"
        );


    // Constructors

        //- Construct from mesh, coordinates and topology
        //  Other properties initialised as null
        inline EvapDropletParcel
        (
            const polyMesh& mesh,
            const barycentric& coordinates,
            const label celli,
            const label tetFacei,
            const label tetPti
        );

        //- Construct from a position and a cell, searching for the rest of the
        //  required topology. Other properties are initialised as null.
        inline EvapDropletParcel
        (
            const polyMesh& mesh,
            const vector& position,
            const label celli
        );

        //- Construct from Istream
        EvapDropletParcel
        (
            const polyMesh& mesh,
            Istream& is,
            bool readFields = true,
            bool newFormat = true
        );

        //- Construct as a copy
        EvapDropletParcel
        (
            const EvapDropletParcel& p,
            const polyMesh& mesh
        );

        //- Construct as a copy
        EvapDropletParcel(const EvapDropletParcel& p);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone() const
        {
            return autoPtr<particle>(new EvapDropletParcel<ParcelType>(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<particle> clone(const polyMesh& mesh) const
        {
            return autoPtr<particle>
            (
                new EvapDropletParcel<ParcelType>(*this, mesh)
            );
        }

        //- Factory class to read-construct particles used for
        //  parallel transfer
        class iNew
        {
            const polyMesh& mesh_;

        public:

            iNew(const polyMesh& mesh)
            :
                mesh_(mesh)
            {}

            autoPtr<EvapDropletParcel<ParcelType>> operator()(Istream& is) const
            {
                return autoPtr<EvapDropletParcel<ParcelType>>
                (
                    new EvapDropletParcel<ParcelType>(mesh_, is, true)
                );
            }
        };


    // Member Functions

        // Access

            //- Return const access to initial mass [kg]
            inline scalar mass0() const;

            inline scalar mDot() const;

        // Edit

            //- Return access to initial mass [kg]
            inline scalar& mass0();

        // Main calculation loop

            //- Set cell values
            template<class TrackCloudType>
            void setCellValues(TrackCloudType& cloud, trackingData& td);

            //- Correct cell values using latest transfer information
            template<class TrackCloudType>
            void cellValueSourceCorrection
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );

            //- Update parcel properties over the time interval
            template<class TrackCloudType>
            void calc
            (
                TrackCloudType& cloud,
                trackingData& td,
                const scalar dt
            );


        // I-O

            //- Read - no composition
            template<class CloudType>
            static void readFields(CloudType& c);

            //- Write - no composition
            template<class CloudType>
            static void writeFields(const CloudType& c);

            //- Write individual parcel properties to stream
            void writeProperties
            (
                Ostream& os,
                const wordRes& filters,
                const word& delim,
                const bool namesOnly = false
            ) const;

            //- Read particle fields as objects from the obr registry
            //  - no composition
            template<class CloudType>
            static void readObjects
            (
                CloudType& c,
                const objectRegistry& obr
            );

            //- Write particle fields as objects into the obr registry
            //  - no composition
            template<class CloudType>
            static void writeObjects
            (
                const CloudType& c,
                objectRegistry& obr
            );

    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const EvapDropletParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "EvapDropletParcelI.H"
#include "EvapDropletParcelTrackingDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EvapDropletParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
