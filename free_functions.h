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
#include "Archipelago.h"

double getLogGrowth(const Island& island) noexcept;

double extractSumOfRates(const Island& island);


#endif //ARCHIDAISIE_FREE_FUNCTIONS_H
