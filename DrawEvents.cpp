//
// Created by Bastophiles on 06.09.2018.
//

#include "DrawEvents.h"

// --- distribution functions ---
int drawDisEvent(const std::vector<double> &vecRates, std::mt19937_64 &prng) {
//draw a discrete distribution event
    std::discrete_distribution<int> drawEvent(vecRates.begin(), vecRates.end());
    return drawEvent(prng);
}

int drawUniEvent(const int &botBoundary, const int &topBoundary, std::mt19937_64 &prng) {
//draw a uniform distribution event
    // using "topBoundary" (= size of vector - 1 = position of last element of vector)
    // because this function works with different types of vectors: vector<int> & vector<pair<int,int> >
    std::uniform_int_distribution<int> drawEvent(botBoundary, topBoundary);
    return drawEvent(prng);
}
