/*---------------------------------------------------------------------------*\
                                       8888888888                              
                                       888                                     
                                       888                                     
  88888b.d88b.  88888b.d88b.   .d8888b 8888888  .d88b.   8888b.  88888b.d88b.  
  888 "888 "88b 888 "888 "88b d88P"    888     d88""88b     "88b 888 "888 "88b 
  888  888  888 888  888  888 888      888     888  888 .d888888 888  888  888 
  888  888  888 888  888  888 Y88b.    888     Y88..88P 888  888 888  888  888 
  888  888  888 888  888  888  "Y8888P 888      "Y88P"  "Y888888 888  888  888 
------------------------------------------------------------------------------- 

License
    This file is part of mmcFoam.

    mmcFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    mmcFoam is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with mmcFoam. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "particleThermoMixture.H"
#include "fvMesh.H"


template<class ThermoType>
Foam::List<Foam::word>
Foam::particleThermoMixture<ThermoType>::constructSpecies
(
    const dictionary& thermoDict
)
{
    List<word> species;
    if 
    (
        thermoDict.subDict("thermoType").get<word>("mixture") 
     != "reactingMixture"
    )
    {
        return thermoDict.get<wordList>("species");
    }
    return species;
}


template<class ThermoType>
const ThermoType& Foam::particleThermoMixture<ThermoType>::constructSpeciesThermo
(
    const dictionary& thermoDict
)
{
    speciesData_.resize(species_.size());

    if 
    (
        thermoDict.subDict("thermoType").get<word>("mixture") 
     == "reactingMixture"
    )
    {
        const ReactionTable<ThermoType>& thermoData = reader_->speciesThermo();
        forAll(species_, i)
        {
            speciesData_.set
            (
                i,
                new ThermoType(*thermoData[species_[i]])
            );
        }
    }
    else
    {
        forAll(species_, i)
        {
            speciesData_.set
            (
                i,
                new ThermoType(thermoDict.subDict(species_[i]))
            );
        }


    }

    return speciesData_[0];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::particleThermoMixture<ThermoType>::particleThermoMixture
(
    const dictionary& thermoDict
)
:
    baseParticleThermoMixture(thermoDict),
    species_(constructSpecies(thermoDict)),
    reader_
    (
        thermoDict.subDict("thermoType").get<word>("mixture")
     == "reactingMixture" ? 
        chemistryReader<ThermoType>::New(thermoDict, species_)
        : 
        nullptr
    ),
    mixture_("mixture",constructSpeciesThermo(thermoDict))
{
    reader_.clear();
}


template<class ThermoType>
Foam::particleThermoMixture<ThermoType>::particleThermoMixture
(
    const particleThermoMixture<ThermoType>& pModel
)
:
    baseParticleThermoMixture(pModel),
    species_(pModel.species_),
    reader_(),
    speciesData_(pModel.speciesData_),
    mixture_(pModel.mixture_)
{}



template<class ThermoType>
Foam::autoPtr<Foam::baseParticleThermoMixture>
Foam::particleThermoMixture<ThermoType>::clone() const
{
    return autoPtr<baseParticleThermoMixture>
    (
        new particleThermoMixture(*this)
    );
}


template<class ThermoType>
void Foam::particleThermoMixture<ThermoType>::update
(
    const List<scalar>& Y
) const
{
    mixture_ = Y[0]*speciesData_[0];
  
    for (label n=1; n<Y.size(); n++)
    {
        mixture_ += Y[n]*speciesData_[n];
    }
}


template<class ThermoType>
Foam::scalar Foam::particleThermoMixture<ThermoType>::mu
(
    const scalar p, const scalar T
) const
{
    return mixture_.mu(p,T);
}


template<class ThermoType>
Foam::scalar Foam::particleThermoMixture<ThermoType>::rho
(
    const scalar p, const scalar T
) const
{
    return mixture_.rho(p,T);
}


template<class ThermoType>
Foam::scalar Foam::particleThermoMixture<ThermoType>::Hs
(
    const scalar p, const scalar T
) const
{
    return mixture_.Hs(p,T);
}


template<class ThermoType>
Foam::scalar Foam::particleThermoMixture<ThermoType>::Cp
(
    const scalar p, const scalar T
) const
{
    return mixture_.Cp(p,T);
}


template<class ThermoType>
Foam::scalar Foam::particleThermoMixture<ThermoType>::kappa
(
    const scalar p, const scalar T
) const
{
    return mixture_.kappa(p,T);
}

// ************************************************************************* //
