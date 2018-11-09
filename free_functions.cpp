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

Species findOldestSpecies(const std::vector<Species> speciesVec)
{
    Species oldest = Species();
    for (auto& species : speciesVec) {
        if (species.getBirth() > oldest.getBirth())
            oldest = species;
    }
    return oldest;
}

STTtable mergeSTTtables(const vector<STTtable>& STTvec, const int n_samples)
{
    double oldest = STTvec[0].getSTTtable()[0].getTime();
    for (auto& stt : STTvec) {
        assert(stt.getSTTtable()[0].getTime()==oldest);
    }
    STTtable fullSTT = STTtable(n_samples);
    return fullSTT;
}

