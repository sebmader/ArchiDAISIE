//
// Created by Bastophiles on 08.11.2018.
//

#include "STTtable.h"

using namespace std;

STTtable::STTtable(const STT& firstRow)
{
    mSTTtable = vector<STT>(1,firstRow);
}

void STTtable::updateSTTtable(const Archipelago& archi, const double& time)
{
    int nImmigrants = 0;
    int nAnagenetic = 0;
    int nCladogenetic = 0;
    vector<Island> islands = archi.getIslands();
    for (auto& isl : islands) {
        vector<Species> species = isl.getSpecies();
        for (auto& sp : species) {
            switch(sp.getStatus()) {
            case 'I':
                ++nImmigrants;
                break;
            case 'A':
                ++nAnagenetic;
                break;
            case 'C':
                ++nCladogenetic;
                break;
            default:
                assert(!"Shouldn't get here!\n"); //!OCLINT
                throw logic_error("Status of species is unknown.\n");
            }
        }
    }
    mSTTtable.emplace_back(time,nImmigrants,nAnagenetic,nCladogenetic);
}

unsigned long STTtable::size()
{
    return mSTTtable.size();
}

const vector<STT>& STTtable::getSTTtable() const
{
    return mSTTtable;
}

