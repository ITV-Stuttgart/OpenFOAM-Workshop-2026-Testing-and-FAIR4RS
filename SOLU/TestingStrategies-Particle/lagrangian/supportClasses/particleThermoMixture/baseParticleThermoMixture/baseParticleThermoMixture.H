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
    Foam::baseParticleThermoMixture

Description
    Abstract base class to allow the computation of thermophysical properties
    based on a given species mass fraction composition, pressure and 
    temperature


SourceFiles
    baseParticleThermoMixture.C

\*---------------------------------------------------------------------------*/

#ifndef baseParticleThermoMixture_H
#define baseParticleThermoMixture_H

#include "dictionary.H"
#include "scalarList.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class baseParticleMixture Declaration
\*---------------------------------------------------------------------------*/

class baseParticleThermoMixture
{
public:
    
    TypeName("baseParticleThermoMixture");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        baseParticleThermoMixture,
        dictionary,
        (const dictionary& thermoDict),
        (thermoDict)
    );

    //- Selector
    static autoPtr<baseParticleThermoMixture> New(const dictionary& thermoDict);
    

    // Constructors

        //- Construct from dictionary, mesh and phase name
        baseParticleThermoMixture(const dictionary&);

        //- Construct as copy
        baseParticleThermoMixture(const baseParticleThermoMixture&);

        virtual ~baseParticleThermoMixture() {}

    // Clone

        virtual autoPtr<baseParticleThermoMixture> clone() const;

    // Member functions

        //- Update the mixture composition based on the mass fractions
        virtual void update(const List<scalar>& Y) const = 0;


        // Calculate thermodynamic states mixture

            //- Viscosity
            inline virtual scalar mu(const scalar p, const scalar T) const = 0;

            //- Density 
            inline virtual scalar rho(const scalar p, const scalar T) const = 0;

            //- Enthalpy
            inline virtual scalar Hs(const scalar p, const scalar T) const = 0;

            //- Heat capacity
            inline virtual scalar Cp(const scalar p, const scalar T) const = 0;

            //- Heat conductivity
            inline virtual scalar kappa(const scalar p, const scalar T) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
