//
// Created by Bastophiles on 29.09.2018.
//

#include "free_functions.h"

using namespace std;

double getLogGrowth(const Island& island)
{
    return 1.0 - static_cast<double>(island.getNSpecies())
            / island.getCarryingCap();
}

double extractSumOfRates(const Island& island) noexcept
{
    vector<double> localRates = island.getLocalRates();
    double sumRates = 0;
    for (double rate : localRates)
        sumRates += rate;
    return sumRates;
}

