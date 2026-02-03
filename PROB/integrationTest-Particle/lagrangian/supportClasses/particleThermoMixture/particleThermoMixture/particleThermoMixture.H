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

Class
    Foam::particleThermoMixture

Description
    Template class for the particle thermo model.
    Constructs the thermophysical model based on the given thermo dict. 

SourceFiles
    particleThermoMixture.C

\*---------------------------------------------------------------------------*/

#ifndef particleThermoMixture_H
#define particleThermoMixture_H

#include "baseParticleThermoMixture.H"
#include "speciesTable.H"
#include "chemistryReader.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class particleThermoMixture Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class particleThermoMixture
:
    public baseParticleThermoMixture
{
public:
    //- The type of thermodynamic
    typedef ThermoType thermoType;

private:
    // Private data
        
        //- Table of specie names
        speciesTable species_;

        //- The chemistry reader
        autoPtr<chemistryReader<ThermoType>> reader_;

        PtrList<thermoType> speciesData_;

        mutable thermoType mixture_;

        //- Create thermoType models for the species
        List<word> constructSpecies(const dictionary& dict);

        //- Create thermoType models for the species
        const ThermoType& constructSpeciesThermo(const dictionary& dict);



public:

    // Constructors

        //- Construct from dictionary
        particleThermoMixture(const dictionary&);

        //- Construct as copy
        particleThermoMixture(const particleThermoMixture<ThermoType>&);

    // Member functions

        TypeName("particleThermoMixture");

        const thermoType& mixture() const
        {
            return mixture_;
        }

        virtual autoPtr<baseParticleThermoMixture> clone() const override;


        //- Update the mixture composition based on the mass fractions
        void update(const List<scalar>& Y) const;

        // Calculate thermodynamic states mixture

            //- Viscosity
            inline scalar mu(const scalar p, const scalar T) const;

            //- Density 
            inline scalar rho(const scalar p, const scalar T) const;

            //- Enthalpy
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Heat capacity
            inline virtual scalar Cp(const scalar p, const scalar T) const;

            //- Heat conductivity
            inline virtual scalar kappa(const scalar p, const scalar T) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "particleThermoMixture.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
