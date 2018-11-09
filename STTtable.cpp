//
// Created by Sebastian Mader on 08.11.2018.
//

#include "STTtable.h"

using namespace std;

STTtable::STTtable(const int size)
{
    assert(size >= 0);
    mSTTtable = vector<STT>(static_cast<unsigned>(size));
}

void STTtable::updateSTTtable(const Archipelago& archi, const double& time)
{
    int nImmigrants = 0;
    int nAnagenetic = 0;
    int nCladogenetic = 0;
    vector<Species> lineages;
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
            bool lineageIsPresent = false;
            for (auto& lineage : lineages) {
                if (lineage.isSister(sp))
                    lineageIsPresent = true;
            }
            if (!lineageIsPresent)
                lineages.push_back(sp);
        }
    }
    int nColonisations = static_cast<int>(lineages.size());
    mSTTtable.emplace_back(time,nImmigrants,nAnagenetic,nCladogenetic,nColonisations);
}

unsigned long STTtable::size()
{
    return mSTTtable.size();
}

const vector<STT>& STTtable::getSTTtable() const noexcept
{
    return mSTTtable;
}

ostream& operator<<(ostream& os, const STTtable& table)
{
    os << "time" << ',' << "nI" << ',' << "nA" << ',' << "nC" << '\n';
    for (auto& stt : table.getSTTtable()) {
        os << stt.getTime() << ',' << stt.getNImmigrants() << ','
           << stt.getNAnagenetic() << ',' << stt.getNCladogenetic() << '\n';
    }
    os << '\n';
    return os;
}

