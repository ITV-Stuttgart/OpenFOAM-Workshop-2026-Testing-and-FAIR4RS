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

Description
    Test the pairing method of the mixParticleModel

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2022
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "Cloud.H"
#include "mixParticleModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Create a dummy class to expose the buildParticleList and findPairs
// function of the mixParticleModel
namespace Foam
{
/*---------------------------------------------------------------------------*\
                Class dummyMixModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class dummyMixModel
:
    public mixParticleModel<CloudType>
{
    // Typedef for particle type
    typedef typename CloudType::particleType particleType;
    typedef typename mixParticleModel<CloudType>::eulerianFieldData eulerianFieldData;
    
    // Private data
    
    public:

        //- Runtime type information
        TypeName("dummyMixModel");

    // Constructors

        //- Construct from dictionary
        dummyMixModel(const dictionary& dict, CloudType& owner, const mmcVarSet& Xi)
        : mixParticleModel<CloudType>(dict,owner, typeName, Xi)
        {}

        //- Construct copy
        dummyMixModel(const dummyMixModel<CloudType>& cm) 
        : mixParticleModel<CloudType>(cm) {}

        //- Construct and return a clone
        virtual autoPtr<CloudMixingModel<CloudType> > clone() const
        {
            return autoPtr<CloudMixingModel<CloudType> >
            (
                new dummyMixModel<CloudType>(*this)
            );
        }


    // Destructor
        virtual ~dummyMixModel() = default;
        

    // Member Functions
    
        //- Dummy function
        void mixpair
        (
            particleType& p,
            const eulerianFieldData& pEulerianFields,
            particleType& q,
            const eulerianFieldData& qEulerianFields,
            scalar& deltaT
        ) {};

        //- Returns the pairing as particle original IDs
        void callFindPairs
        (
            DynamicList<List<label>>& particlePairs
        )
        {
            this->buildParticleList(-GREAT,GREAT);

            DynamicList<List<label>> particlePairsTemp = this->particlePairs_;

            for (const auto& pair : particlePairsTemp)
            {
                List<label> temp = pair;
                bool isLocal = false;
                forAll(pair,j)
                {
                    if (this->eulerianFieldDataList_[pair[j]].local())
                        isLocal = true;
                }
                if (isLocal)
                {
                    forAll(pair,j)
                    {
                        temp[j] = this->particleList_
                        [
                            this->eulerianFieldDataList_[pair[j]].particleIndex()
                        ].wt();
                    }
                
                    particlePairs.append(temp);
                }
            }
        }


};

} // End of namespace Foam

