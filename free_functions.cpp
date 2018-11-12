//
// Created by Sebastian Mader on 29.09.2018.
//

#include "free_functions.h"

using namespace std;

double getLogGrowth(const Island& island) noexcept
{
    return 1.0 - static_cast<double>(island.getNSpecies())
            / static_cast<double>(island.getCarryingCap());
}

double extractSumOfRates(const Island& island)
{
    vector<double> localRates = island.getLocalRates();
    double sumRates = 0;
    for (double rate : localRates)
        sumRates += rate;
    return sumRates;
}

Species findOldestSpecies(const vector<Species>& speciesVec)
{
    Species oldest = Species();
    for (auto& species : speciesVec) {
        if (species.getBirth() > oldest.getBirth())
            oldest = species;
    }
    return oldest;
}

STTtable mergeSTTtables(const std::vector<STTtable>& STTvec, const int& n_samples)
{
    double oldestAge = STTvec[0].getSTTtable()[0].getTime();
    for (auto& stt : STTvec) {
        assert(stt.getSTTtable()[0].getTime()==oldestAge);
    }
    STTtable fullSTT = STTtable((unsigned) n_samples+1, STT(oldestAge, 0, 0, 0, 0));

    // count species states per predefined time slice (island age / sample freq)
    for(size_t i = 1; i < fullSTT.size(); ++i) {
        double time = oldestAge - (oldestAge/n_samples) * i;
        assert(time >= 0.0);
        int nImmi = 0;
        int nAna = 0;
        int nClado = 0;
        int nColo = 0;
        for (auto& sttTab : STTvec) {
            for (int j = (int)sttTab.size()-1; j >= 0; --j) {
                if (sttTab[j].getTime() >= time) {
                    nImmi += sttTab[j].getNImmigrants();
                    nAna += sttTab[j].getNAnagenetic();
                    nClado += sttTab[j].getNCladogenetic();
                    nColo += sttTab[j].getNColonisations();
                    break;
                }
            }
        }
        fullSTT[i] = STT(time,nImmi,nAna,nClado,nColo);
    }
    return fullSTT;
}

