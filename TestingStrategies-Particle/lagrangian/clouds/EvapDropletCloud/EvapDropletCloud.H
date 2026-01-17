/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::EvapDropletCloud

Description
    Templated base class for flashing droplet cloud

    - Adds to thermo cloud 

SourceFiles
    EvapDropletCloudI.H
    EvapDropletCloud.C

\*---------------------------------------------------------------------------*/

#ifndef EvapDropletCloud_H
#define EvapDropletCloud_H

#include "ThermoCloud.H"
#include "evapDropletCloud.H"
#include "baseParticleThermoMixture.H"

// Submodels
#include "EvapModel.H"
#include "liquidProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class EvapModel;

/*---------------------------------------------------------------------------*\
                        Class EvapDropletCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class EvapDropletCloud
:
    public CloudType,
    public evapDropletCloud
{
public:

    // Public typedefs

        //- Type of cloud this cloud was instantiated for
        typedef CloudType cloudType;

        //- Type of parcel the cloud was instantiated for
        typedef typename CloudType::particleType parcelType;

        //- Convenience typedef for this cloud type
        typedef EvapDropletCloud<CloudType> EvapDropletCloudType;


private:

    // Private data

        //- Cloud copy pointer
        autoPtr<EvapDropletCloud<CloudType>> cloudCopyPtr_;


    // Private member functions

        //- Disallow default bitwise copy construct
        EvapDropletCloud(const EvapDropletCloud&);

        //- Disallow default bitwise assignment
        void operator=(const EvapDropletCloud&);


protected:

    // Protected data

        //- SLG thermodynamics package
        const SLGThermo& thermo_;

        //- Parcel constant properties
        typename parcelType::constantProperties constProps_;

        //- Name of the droplet species -- only single component
        word speciesName_;

        //- Index of the species in the gas composition
        label speciesIndex_;

        // References to the cloud sub-models
        
            //-  Evaporation model
            autoPtr<EvapModel<EvapDropletCloudType>>
                evapModel_;

            //- Liquid properties model 
            autoPtr<liquidProperties> liquidProperties_;

            //- Access to the thermo properties for a user specified mixture
            //  Note: This allows access to the functions without modifying the 
            //        state of the mixture. Otherwise the carrier() function 
            //        of SLGThermo can be used. 
            autoPtr<baseParticleThermoMixture> thermoMixture_;

        // Sources
 
            //- Mass transfer fields - one per carrier phase specie
            autoPtr<volScalarField::Internal> rhoTrans_;

    // Protected Member Functions

         // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Reset state of cloud
            void cloudReset(EvapDropletCloud<CloudType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        EvapDropletCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            const SLGThermo& thermo,
            bool readFields = true
        );

        //- Copy constructor with new name
        EvapDropletCloud(EvapDropletCloud<CloudType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        EvapDropletCloud
        (
            const fvMesh& mesh,
            const word& name,
            const EvapDropletCloud<CloudType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> clone(const word& name)
        {
            return autoPtr<Cloud<parcelType>>
            (
                new EvapDropletCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
        {
            return autoPtr<Cloud<parcelType>>
            (
                new EvapDropletCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~EvapDropletCloud();

    // Member Functions
 
        // Access
 
            //- Return a reference to the cloud copy
            inline const EvapDropletCloud& cloudCopy() const;
 
            //- Return the constant properties
            inline const typename parcelType::constantProperties&
                constProps() const;
 
            //- Return access to the constant properties
            inline typename parcelType::constantProperties& constProps();

            //- Return the species name 
            inline const word& dropletSpeciesName() const;

            //- Return the species index
            inline const label& dropletSpeciesIndex() const;

            //- Return the SLG thermo model
            inline const SLGThermo& thermo() const;

            //- Return the particle mixture state
            inline const baseParticleThermoMixture& thermoMixture() const;

            // Sub-models
 
                //- Return const access to reacting composition model
                inline const EvapModel<EvapDropletCloudType>&
                    evapModel() const;
 
                //- Return reference to reacting phase change model
                inline EvapModel<EvapDropletCloudType>&
                    evapModel();

                //- Additional fluid properties, e.g., sat pressure, 
                //  diffusion coefficient, latent heat of vaporization etc. 
                inline const liquidProperties& liquid() const;

                //- Additional fluid properties, e.g., sat pressure, 
                //  diffusion coefficient, latent heat of vaporization etc. 
                inline liquidProperties& liquid();
 
 
        // Sources

            //- Mass

                //- Return reference to mass source
                inline volScalarField::Internal&
                    rhoTrans();

                //- Return const reference to mass source
                inline const volScalarField::Internal&
                    rhoTrans() const;

                //- Return tmp total mass source for carrier phase
                //  - fully explicit
                inline tmp<volScalarField::Internal> Srho() const;

                //- Return total mass source term [kg/m3/s]
                inline tmp<fvScalarMatrix> Srho(volScalarField& rho) const;

                //- Return mass source term for specie i - specie eqn
                inline tmp<fvScalarMatrix> SYi
                (
                    const label i,
                    volScalarField& Yi
                ) const;


    // Cloud evolution functions

        //- Set parcel thermo properties
        void setParcelThermoProperties
        (
            parcelType& parcel,
            const scalar lagrangianDt
        );

        //- Check parcel properties
        void checkParcelProperties
        (
            parcelType& parcel,
            const scalar lagrangianDt,
            const bool fullyDescribed
        );

        //- Store the current cloud state
        void storeState();

        //- Reset the current cloud to the previously stored state
        void restoreState();

        //- Reset the cloud source terms
        void resetSourceTerms();

        //- Apply relaxation to (steady state) cloud sources
        void relaxSources(const EvapDropletCloud<CloudType>& cloudOldTime);

        //- Apply scaling to (transient) cloud sources
        void scaleSources();

        //- Evolve the cloud
        void evolve();


    // Mapping

        //- Remap the cells of particles corresponding to the
        //  mesh topology change with a default tracking data object
        virtual void autoMap(const mapPolyMesh&);

    // I-O

        //- Print cloud information
        void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "EvapDropletCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "EvapDropletCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //