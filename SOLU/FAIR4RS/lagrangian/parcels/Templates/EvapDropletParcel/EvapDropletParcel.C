/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::EvapDropletParcel<ParcelType>::EvapDropletParcel
(
    const EvapDropletParcel<ParcelType>& p
)
:
    ParcelType(p),
    mass0_(p.mass0_)
{}


template<class ParcelType>
Foam::EvapDropletParcel<ParcelType>::EvapDropletParcel
(
    const EvapDropletParcel<ParcelType>& p,
    const polyMesh& mesh
)
:
    ParcelType(p, mesh),
    mass0_(p.mass0_)
{}


// * * * * * * * * * * * * *  Member Functions * * * * * * * * * * * * * * * //


template<class ParcelType>
template<class TrackCloudType>
void Foam::EvapDropletParcel<ParcelType>::setCellValues
(
    TrackCloudType& cloud,
    trackingData& td
)
{
    ParcelType::setCellValues(cloud, td);

    td.pc() = td.pInterp().interpolate
    (
        this->coordinates(),
        this->currentTetIndices()
    );

    List<scalar>& Y = td.Yc();

    forAll(Y,i)
    {
        Y[i] = td.YInterp(i).interpolate
        (
            this->coordinates(),
            this->currentTetIndices()
        );
    }
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::EvapDropletParcel<ParcelType>::cellValueSourceCorrection
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    scalar maxMassI = 0.0;

    scalar dm = cloud.rhoTrans()[this->cell()];
    maxMassI = max(maxMassI, mag(dm));
    scalar addedMass = dm;

    if (maxMassI < ROOTVSMALL)
    {
        return;
    }

    const scalar massCell = this->massCell(td);

    td.rhoc() += addedMass/cloud.pMesh().cellVolumes()[this->cell()];

    const scalar massCellNew = massCell + addedMass;
    td.Uc() = (td.Uc()*massCell + cloud.UTrans()[this->cell()])/massCellNew;
}


template<class ParcelType>
template<class TrackCloudType>
void Foam::EvapDropletParcel<ParcelType>::calc
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar dt
)
{
    // Reference Conditions
    // ~~~~~~~~~~~~~~~~~~~~

    // Define local properties at beginning of time step
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 
    const scalar np0 = this->nParticle_;
    const scalar d0 = this->d_;
    const vector& U0 = this->U_;
    const scalar T0 = this->T_;
    const scalar mass0 = this->mass();


    // Calc surface values based on 2/3 rule -- required for fluid properties
    scalar Ts, rhos, mus, Prs, kappas;
    this->calcSurfaceValues(cloud, td, T0, Ts, rhos, mus, Prs, kappas);

    // Calculate characteristic numbers for the surface properties
    const scalar Res = this->Re(rhos, U0, td.Uc(), d0, mus);

    // Sources
    // ~~~~~~~

    // Explicit momentum source for particle
    vector Su = Zero;

    // Linearised momentum source coefficient
    scalar Spu = 0.0;

    // Momentum transfer from the particle to the carrier phase
    vector dUTrans = Zero;

    // Linearised enthalpy source coefficient
    scalar Sph = 0.0;


    // 1. Compute models that contribute to mass transfer - U, T held constant
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Phase change
    // ~~~~~~~~~~~~

    scalar mDot, QDot;
    cloud.evapModel().mDotAndQDot
    (
        td,
        T0,
        this->U(),
        d0,
        mass0,
        this->Cp_,
        mDot,
        QDot
    );

    // Store for post-processing
    mDot_ = mDot;

    // Calcualte new droplet mass
    scalar mass1 = mass0 - mDot*dt;

    // 2. Update the parcel properties due to change in mass
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    // Remove the particle when mass falls below minimum threshold
    if (np0*mass1 < cloud.constProps().minParcelMass())
    {
        td.keepParticle = false;

        if (cloud.solution().coupled())
        {
            scalar dm = np0*mass0;

            cloud.rhoTrans()[this->cell()] += dm;
            cloud.UTrans()[this->cell()] += dm*U0;
        }

        return;
    }

    // Update the mass
    this->mass0_ = mass1;


    // Calculate the deltaT of droplet temperature
    const scalar deltaT = QDot*(dt/(mass0*this->Cp_));
    this->T_ =std::max(T0 + deltaT,10.0);

    // Update the diameter
    this->d_ = cbrt(mass1/this->rho_*6/pi);

    // Update the heat capacity
    this->Cp_ = cloud.liquid().Cp(td.pc(),this->T_);


    // Motion
    // ~~~~~~

    // Calculate new particle velocity
    this->U_ =
        this->calcVelocity(cloud, td, dt, Res, mus, mass1, Su, dUTrans, Spu);


    // 4. Accumulate carrier phase source terms
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (cloud.solution().coupled())
    {
        scalar dm = np0*(mass0-mass1);

        // Transfer the mass to the Eulerian cell
        cloud.rhoTrans()[this->cell()] += dm;

        // Update momentum transfer
        cloud.UTrans()[this->cell()] += np0*dUTrans;
        cloud.UCoeff()[this->cell()] += np0*Spu;

        // Vapor is generated at saturation temperature
        const scalar TSat = cloud.liquid().pvInvert(td.pc());
        const scalar hs = 
            cloud.thermo().carrier().Hs(cloud.dropletSpeciesIndex(),td.pc(),TSat);

        // Update sensible enthalpy transfer
        cloud.hsTrans()[this->cell()] += dm*hs;
        cloud.hsCoeff()[this->cell()] += np0*Sph;
    }
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "EvapDropletParcelIO.C"

// ************************************************************************* //