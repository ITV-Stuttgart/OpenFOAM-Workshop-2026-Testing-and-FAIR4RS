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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::EvapDropletCloud<CloudType>&
Foam::EvapDropletCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}


template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::EvapDropletCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::EvapDropletCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
const Foam::word& 
Foam::EvapDropletCloud<CloudType>::dropletSpeciesName() const
{
    return speciesName_;
}


template<class CloudType>
const Foam::label& 
Foam::EvapDropletCloud<CloudType>::dropletSpeciesIndex() const
{
    return speciesIndex_;
}


template<class CloudType>
const Foam::SLGThermo& 
Foam::EvapDropletCloud<CloudType>::thermo() const
{
    return thermo_;
}


template<class CloudType>
const Foam::baseParticleThermoMixture& 
Foam::EvapDropletCloud<CloudType>::thermoMixture() const
{
    return thermoMixture_;
}


template<class CloudType>
inline const Foam::EvapModel<Foam::EvapDropletCloud<CloudType>>&
Foam::EvapDropletCloud<CloudType>::evapModel() const
{
    return *evapModel_;
}


template<class CloudType>
inline Foam::EvapModel<Foam::EvapDropletCloud<CloudType>>&
Foam::EvapDropletCloud<CloudType>::evapModel()
{
    return *evapModel_;
}


template<class CloudType>
inline const Foam::liquidProperties&
Foam::EvapDropletCloud<CloudType>::liquid() const
{
    return *liquidProperties_;
}


template<class CloudType>
inline Foam::liquidProperties&
Foam::EvapDropletCloud<CloudType>::liquid()
{
    return *liquidProperties_;
}

template<class CloudType>
inline Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::EvapDropletCloud<CloudType>::rhoTrans()
{
    return rhoTrans_.ref();
}


template<class CloudType>
inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::EvapDropletCloud<CloudType>::rhoTrans() const
{
    return rhoTrans_();
}


template<class CloudType>
inline Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>
Foam::EvapDropletCloud<CloudType>::Srho() const
{
    tmp<volScalarField::Internal> tRhoi
    (
        new volScalarField::Internal
        (
            IOobject
            (
                this->name() + ":rhoTrans",
                this->db().time().timeName(),
                this->db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            this->mesh(),
            dimensionedScalar
            (
                rhoTrans_().dimensions()/dimTime/dimVolume, Zero
            )
        )
    );

    if (this->solution().coupled())
    {
        scalarField& rhoi = tRhoi.ref();
        rhoi = rhoTrans_()/(this->db().time().deltaTValue()*this->mesh().V());
    }

    return tRhoi;
}


template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix>
Foam::EvapDropletCloud<CloudType>::Srho(volScalarField& rho) const
{
    if (this->solution().coupled())
    {
        tmp<volScalarField> trhoTrans
        (
            new volScalarField
            (
                IOobject
                (
                    this->name() + ":rhoTrans",
                    this->db().time().timeName(),
                    this->db(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                this->mesh(),
                dimensionedScalar(dimMass/dimTime/dimVolume, Zero)
            )
        );

        scalarField& sourceField = trhoTrans.ref();

        if (this->solution().semiImplicit("rho"))
        {

            sourceField += rhoTrans_();
            sourceField /= this->db().time().deltaTValue()*this->mesh().V();

            return fvm::SuSp(trhoTrans()/rho, rho);
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(rho, dimMass/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();

            sourceField += rhoTrans_();

            fvm.source() = -trhoTrans()/this->db().time().deltaT();

            return tfvm;
        }
    }

    return tmp<fvScalarMatrix>::New(rho, dimMass/dimTime);
}

template<class CloudType>
inline Foam::tmp<Foam::fvScalarMatrix> 
Foam::EvapDropletCloud<CloudType>::SYi
(
    const label i,
    volScalarField& Yi
) const
{
    if (this->solution().coupled() && i == speciesIndex_)
    {
        if (this->solution().semiImplicit("Yi"))
        {
            tmp<volScalarField> trhoTrans
            (
                new volScalarField
                (
                    IOobject
                    (
                        this->name() + ":rhoTrans",
                        this->db().time().timeName(),
                        this->db(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    this->mesh(),
                    dimensionedScalar(dimMass/dimTime/dimVolume, Zero)
                )
            );
 
            volScalarField& sourceField = trhoTrans.ref();
 
            sourceField.primitiveFieldRef() =
                rhoTrans_()/(this->db().time().deltaTValue()*this->mesh().V());
 
            const dimensionedScalar YiSMALL("YiSMALL", dimless, SMALL);
 
            return
                fvm::Sp(neg(sourceField)*sourceField/(Yi + YiSMALL), Yi)
              + pos0(sourceField)*sourceField;
        }
        else
        {
            tmp<fvScalarMatrix> tfvm(new fvScalarMatrix(Yi, dimMass/dimTime));
            fvScalarMatrix& fvm = tfvm.ref();
 
            fvm.source() = -rhoTrans_()/this->db().time().deltaTValue();
 
            return tfvm;
        }
    }
 
    return tmp<fvScalarMatrix>::New(Yi, dimMass/dimTime);
}
// ************************************************************************* //