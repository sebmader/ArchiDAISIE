//
// Created by Bastophiles on 06.09.2018.
//

#include "DrawEvents.h"

// --- distribution functions ---
int drawDisEvent(const std::vector<double> &vecRates, std::mt19937_64 &prng)
{  //draw a discrete distribution event
    std::discrete_distribution<int> drawEvent(vecRates.begin(), vecRates.end());
    int event = drawEvent(prng);
    return event;
}

int drawUniEvent(const int &botBoundary, const int &topBoundary, std::mt19937_64 &prng)
{  // draw a uniform distribution event
    std::uniform_int_distribution<int> drawEvent(botBoundary, topBoundary);
    int event = drawEvent(prng);
    return event;
}

SpeciesID drawUniEvent(const std::vector<SpeciesID> &speciesVec, std::mt19937_64 &prng)
{  // draw species from vector
    const int vecSize = static_cast<int>(speciesVec.size());
    std::uniform_int_distribution<int> drawSpecies(0, vecSize-1);
    int speciesPos = drawSpecies(prng);
    const SpeciesID specID = speciesVec[speciesPos];
    return specID;
}
