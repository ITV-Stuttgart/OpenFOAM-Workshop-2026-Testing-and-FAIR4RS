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

#include "AbramzonSirignanoModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::AbramzonSirignanoModel<CloudType>::AbramzonSirignanoModel
(
    const dictionary& dict,
    CloudType& owner
)
:
    EvapModel<CloudType>(owner)
{}


template<class CloudType>
Foam::AbramzonSirignanoModel<CloudType>::AbramzonSirignanoModel
(
    const AbramzonSirignanoModel<CloudType>& flashingModel
)
:
    EvapModel<CloudType>(flashingModel.owner_)
{}


template<class CloudType>
void Foam::AbramzonSirignanoModel<CloudType>::mDotAndQDot
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
    const auto& cloud = this->owner_;

    // Infinity pressure
    const scalar pInf = td.pc();

    // Infinity temperature
    const scalar TInf = td.Tc();

    // Calculate the mole fraction of the fuel species at the droplet surface
    const scalar pSat = cloud.liquid().pv(td.pc(),Td);

    // Mole fraction of the fuel at the surface
    const scalar XsFuel = min(pSat/pInf,1.0-1E-08);

    // Droplet species index
    const scalar jFuel = cloud.dropletSpeciesIndex();

    // Get the infinity species composition
    const List<scalar>& YInf = td.Yc();

    // Get the molecular weight vector [g/mol]
    List<scalar> MW(YInf.size());
    forAll(MW,i)
    {
        MW[i] = cloud.thermo().carrier().W(i);
    }

    // Calculate the mole fraction of all species at infinity conditions
    const List<scalar> XInf = YtoX(YInf,MW);

    // Check if XInf is fully saturated
    // If already saturated only heat the droplet
    if (XInf[jFuel] >= 1.0)
    {
        mDot = 0;
        const scalar kappaRef = cloud.liquid().kappa(pInf,Td);
        QDot = 2.0*pi*d*kappaRef*(TInf-Td);
        return;
    }

    // Calculate the mole fractions at the droplet surface, except fuel
    List<scalar> Xs(XInf.size());
    
    forAll(Xs,i)
        Xs[i] = XInf[i]*(1.0-XsFuel)/(1.0-XInf[jFuel]);
    Xs[jFuel] = XsFuel;

    // Now calculate again the mass fractions based on the mole fractions
    List<scalar> Ys = XtoY(Xs,MW);

    // Evaluate the average gas properties
    scalar TRef;    // reference temperature
    scalar kappaRef;    // reference heat conductivity
    scalar CpRef;   // heat capacity
    scalar CpFRef;  // heat capacity of the fuel species
    scalar muRef;
    scalar rhoRef;

    calcReferenceValues
    (
        cloud,
        pInf, Td, TInf, Ys, YInf, jFuel,
        TRef, rhoRef, muRef, kappaRef, CpRef, CpFRef
    );

    // Mol weight of the gas
    scalar MW_Mix = 0;
    forAll(MW,i)
    {
        MW_Mix += MW[i]*XInf[i];
    }
    
    const scalar Df = cloud.liquid().D(pInf,TRef,MW_Mix);


    // Calculate the non-dimensional numbers
    const scalar Sc = muRef/(rhoRef*Df);
    const scalar Pr = CpRef*muRef/(kappaRef);
    const scalar Le = Sc/Pr;
    const scalar Re = CloudType::parcelType::Re(rhoRef,Ud,td.Uc(),d,muRef);
    const scalar Sh = 1.0+pow(1.0+Re*Sc,1.0/3.0)*max(1.0,pow(Re,0.077));
    const scalar Nu = 1.0+pow(1.0+Re*Pr,1.0/3.0)*max(1.0,pow(Re,0.077));;


    // Spalding mass transfer number
    const scalar BM  = (Ys[jFuel]-YInf[jFuel])/(1.0-Ys[jFuel]);
    
    const scalar FM  = pow(1.0+BM,0.7)*log(1.0+BM)/BM;
    const scalar ShMod = 2.0+(Sh-2.0)/FM;
    mDot = pi*d*ShMod*rhoRef*Df*log(1.0+BM);
    scalar BT = SMALL;

    if (fabs(BM)>1e-12)
    {
        scalar phi = (CpFRef/CpRef)*(Sh/Nu)*(1.0/Le);  // initial guess
        const scalar BT0 = pow(1.0+BM,phi)-1.0;
        BT = BT0;

        unsigned int iter = 0;
        const unsigned int iterMax = 100;
        const scalar eps = 1.0e-08;
        for (iter=0; iter<iterMax; ++iter)
        {
            const scalar FT = pow(1.0+BT,0.7)*log(1.0+BT)/BT;
            const scalar NuMod = 2.0+(Nu-2.0)/FT;
            phi = (CpFRef/CpRef)*(ShMod/NuMod)*(1.0/Le);
            const scalar BTNew = pow(1.0+BM,phi)-1.0;
            const scalar res = fabs(BTNew-BT);
            BT = BTNew;
            if (res < eps) 
                break;
        }
        if (iter == iterMax-1)
        {
            Pout << "Warning: Iteration of BT did not converge within " 
                 << iter << " iterations!\n";
            BT = BT0;
        }
        
        const scalar Lv = cloud.liquid().hl(pInf,Td);
        QDot  = mDot*(CpFRef*(TInf-Td)/BT-Lv);
    }
    else
    {
        QDot = 2.0*pi*d*kappaRef*(TInf-Td);  // pure heat transfer
    }
}


template<class CloudType>
void Foam::AbramzonSirignanoModel<CloudType>::calcReferenceValues
(
    const CloudType& cloud,
    const scalar pInf,          // Infinity pressure
    const scalar Ts,            // Surface temperature
    const scalar TInf,          // Infinity temperature
    const List<scalar>& Ys,     // Surface species composition
    const List<scalar>& YInf,   // Infinity species composition
    const label jFuel,       // Index of the fuel mass fraction in Y
    scalar& TRef,               // Reference temperature
    scalar& rhoRef,             // Density at reference conditions
    scalar& muRef,              // Viscosity at reference conditions
    scalar& kappaRef,           // heat conductivity at reference condition
    scalar& CpRef,              // heat capacity at reference conditions
    scalar& CpFRef          // heat capacity at reference conditions (pure liquid mass fraction)
) const
{
    // Calculate the reference temperature and mass fraction
    TRef = (TInf + 2.0*Ts)/3.0;
    
    List<scalar> YRef(Ys);
    
    forAll(YRef,i)
    {
        YRef[i] = (YInf[i] + 2.0*Ys[i])/3.0;
    }

    // Determine the thermophysical properties 
    auto& thermoMixture = cloud.thermoMixture();    

    // Update with the current species mass fractions to get the correct
    // mixture properties
    thermoMixture.update(YRef);

    rhoRef = thermoMixture.rho(pInf,TRef);
    muRef = thermoMixture.mu(pInf,TRef);
    kappaRef = thermoMixture.kappa(pInf,TRef);
    CpRef = thermoMixture.Cp(pInf,TRef);

    CpFRef = cloud.thermo().carrier().Cp(jFuel,pInf,TRef);
}


template<class CloudType>
Foam::List<Foam::scalar> 
Foam::AbramzonSirignanoModel<CloudType>::YtoX
(
    const List<scalar>& Y,
    const List<scalar>& MW
) const
{
    // Calculate temporary mean weighted mole mass
    scalar meanYMW = 0;
    forAll(Y,i)
    {
        meanYMW += Y[i]/MW[i];
    }

    List<scalar> Xi(Y.size());

    forAll(Xi,i)
    {
        Xi[i] = (Y[i]/MW[i])/meanYMW;
    }

    return Xi;
}


template<class CloudType>
Foam::List<Foam::scalar> 
Foam::AbramzonSirignanoModel<CloudType>::XtoY
(
    const List<scalar>& X,
    const List<scalar>& MW
) const
{
    // Calculate temporary mean weighted mole mass
    scalar meanXMW = 0;
    forAll(X,i)
    {
        meanXMW += X[i]*MW[i];
    }

    List<scalar> Y(X.size());

    forAll(Y,i)
    {
        Y[i] = X[i]*MW[i]/meanXMW;
    }

    return Y;
}



// ************************************************************************* //
