//
// Created by Bastophiles on 29.09.2018.
//

#ifndef ARCHIDAISIE_FREE_FUNCTIONS_H
#define ARCHIDAISIE_FREE_FUNCTIONS_H

#include <vector>
#include <cassert>
#include "event_type.h"
#include "SpeciesID.h"
#include "Species.h"
#include "Island.h"

double getLogGrowth(const Island& island) noexcept;

double extractSumOfRates(const Island& island);

Species findOldestSpecies(const std::vector<Species>& speciesVec);

#endif //ARCHIDAISIE_FREE_FUNCTIONS_H
