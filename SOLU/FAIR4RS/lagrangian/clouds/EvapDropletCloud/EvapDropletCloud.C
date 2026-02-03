/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "EvapDropletCloud.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::setModels()
{
    // Note: The liquid properties must be created first as the evaporation model
    //       uses them in the constructor
    liquidProperties_.reset
    (
        liquidProperties::New(this->speciesName_)
    );

    evapModel_.reset
    (
        EvapModel<EvapDropletCloudType>::New
        (
            this->subModelProperties(),
            *this
        ).ptr()
    );

    thermoMixture_.reset
    (
        baseParticleThermoMixture::New(thermo_.thermo())
    );
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::cloudReset(EvapDropletCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    liquidProperties_.reset(c.liquidProperties_.ptr());
    evapModel_.reset(c.evapModel_.ptr());
    thermoMixture_.reset(c.thermoMixture_.ptr());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::EvapDropletCloud<CloudType>::EvapDropletCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields
)
:
    CloudType(cloudName, rho, U, g, thermo, false),
    evapDropletCloud(),
    cloudCopyPtr_(nullptr),
    thermo_(thermo),
    speciesName_(this->subModelProperties().template get<word>("dropletComponent")),
    speciesIndex_(this->thermo().carrier().species()[speciesName_]),
    evapModel_(nullptr)
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this);
            this->deleteLostParticles();
        }
    }

    if (speciesIndex_ == -1)
    {
        // Attempt to find the gas name with an N prepended to the name:
        speciesIndex_ = this->thermo().carrier().species()["N"+speciesName_];

        if (speciesIndex_ == -1)
            FatalError << "Species index for " << speciesName_ << " could not be found"
                << " in " << this->thermo().carrier().species()
                << exit(FatalError);
    }
        

    rhoTrans_.reset
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimMass, Zero)
        )
    );

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::EvapDropletCloud<CloudType>::EvapDropletCloud
(
    EvapDropletCloud<CloudType>& c,
    const word& name
)
:
    CloudType(c, name),
    evapDropletCloud(),
    cloudCopyPtr_(nullptr),
    thermo_(c.thermo_),
    speciesName_(c.speciesName_),
    speciesIndex_(c.speciesIndex_),
    evapModel_(c.evapModel_->clone()),
    liquidProperties_(c.liquidProperties_->clone())
{
    rhoTrans_.reset
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            this->mesh(),
            dimensionedScalar(dimMass, Zero)
        )
    );
}


template<class CloudType>
Foam::EvapDropletCloud<CloudType>::EvapDropletCloud
(
    const fvMesh& mesh,
    const word& name,
    const EvapDropletCloud<CloudType>& c
)
:
    CloudType(mesh, name, c),
    evapDropletCloud(),
    cloudCopyPtr_(nullptr),
    thermo_(c.thermo_),
    speciesName_(c.speciesName_),
    speciesIndex_(c.speciesIndex_),
    evapModel_(c.evapModel_->clone()),
    liquidProperties_(c.liquidProperties_->clone())
{}


template<class CloudType>
Foam::EvapDropletCloud<CloudType>::~EvapDropletCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::setParcelThermoProperties
(
    parcelType& parcel,
    const scalar lagrangianDt
)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::checkParcelProperties
(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);

    // derived information - store initial mass
    parcel.mass0() = parcel.mass();
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset
    (
        static_cast<EvapDropletCloud<CloudType>*>
        (
            clone(this->name() + "Copy").ptr()
        )
    );
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
    rhoTrans_.ref().field() = 0.0;
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::relaxSources
(
    const EvapDropletCloud<CloudType>& cloudOldTime
)
{
    CloudType::relaxSources(cloudOldTime);

    typedef volScalarField::Internal dsfType;


    dsfType& rhoT = rhoTrans_();
    const dsfType& rhoT0 = cloudOldTime.rhoTrans();
    this->relax(rhoT, rhoT0, "rho");
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::scaleSources()
{
    CloudType::scaleSources();

    this->scale(rhoTrans_(), "rho");
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::autoMap(const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::EvapDropletCloud<CloudType>::info()
{
    CloudType::info();
}


// ************************************************************************* //