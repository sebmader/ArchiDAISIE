//
// Created by Sebastian Mader on 29.09.2018.
//

#ifndef ARCHIDAISIE_FREE_FUNCTIONS_H
#define ARCHIDAISIE_FREE_FUNCTIONS_H

#include <vector>
#include <cassert>
#include "Island.h"
#include "STTtable.h"

double getLogGrowth(const Island& island) noexcept;

double extractSumOfRates(const Island& island);

Species findOldestSpecies(const std::vector<Species>& speciesVec);

STTtable mergeSTTtables(const std::vector<STTtable>& STTvec, const int& n_samples);

#endif //ARCHIDAISIE_FREE_FUNCTIONS_H
