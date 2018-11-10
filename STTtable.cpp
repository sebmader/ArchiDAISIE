//
// Created by Bastophiles on 08.11.2018.
//

#include "STTtable.h"

using namespace std;

STTtable::STTtable(const unsigned int& size, const STT& rowContent)
{
    mSTTtable = vector<STT>(size, rowContent);
}

void STTtable::push_back(const STT& newRow)
{
    mSTTtable.push_back(newRow);
}

unsigned long STTtable::size() const
{
    return mSTTtable.size();
}

const vector<STT>& STTtable::getSTTtable() const noexcept
{
    return mSTTtable;
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

ostream& operator<<(ostream& os, const STTtable& table)
{
    os << "time" << ',' << "nI" << ',' << "nA" << ',' << "nC" << '\n';
    for (auto& stt : table.getSTTtable()) {
        os << stt.getTime() << ',' << stt.getNImmigrants() << ','
           << stt.getNAnagenetic() << ',' << stt.getNCladogenetic() << '\n';
    }
    return os;
}

