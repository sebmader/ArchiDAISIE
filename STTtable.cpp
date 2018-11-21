//
// Created by Sebastian Mader on 08.11.2018.
//

#include "STTtable.h"
#include "free_functions.h"

using namespace std;

STTtable::STTtable(const unsigned int& size, const STT& rowContent)
{
    mSTTtable = vector<STT>(size, rowContent);
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
                throw logic_error("Status of species is unknown.\n");
            }
            bool lineagePresent = false;
            for (auto& lineage : lineages) {
                if (lineage.isSister(sp))
                    lineagePresent = true;
            }
            if(!lineagePresent)
                lineages.push_back(sp);
        }
    }
    const int nColonisations = static_cast<int>(lineages.size());
    mSTTtable.emplace_back(time,nImmigrants,nAnagenetic,nCladogenetic,nColonisations);
}

ostream& operator<<(ostream& os, const STTtable& table)
{
    os << "Time" << ',' << "nI" << ',' << "nA" << ',' << "nC" << ',' << "nColo" << '\n';
    for (auto& stt : table.getSTTtable()) {
        os << stt.getTime() << ',' << stt.getNImmigrants() << ','
           << stt.getNAnagenetic() << ',' << stt.getNCladogenetic() << ','
           << stt.getNColonisations() << '\n';
    }
    return os;
}

STT& STTtable::operator[](size_t n)
{
    return mSTTtable[n];
}

const STT& STTtable::operator[](size_t n) const
{
    return mSTTtable[n];
}

void STTtable::updateFullSTTtable(const Archipelago& archi, const double& time)
{
    int nImmigrants = 0;
    int nAnagenetic = 0;
    int nCladogenetic = 0;
    vector<Species> lineages;
    Island mergedIsland = archi.makeArchiTo1Island();
    vector<Species> species = mergedIsland.getSpecies();
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
            throw logic_error("Status of species is unknown.\n");
        }
    }
    mSTTtable.emplace_back(time,nImmigrants,nAnagenetic,nCladogenetic,howManyLineages(species));
}