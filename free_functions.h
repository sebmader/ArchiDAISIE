//
// Created by Sebastian Mader on 29.09.2018.
//

#ifndef ARCHIDAISIE_FREE_FUNCTIONS_H
#define ARCHIDAISIE_FREE_FUNCTIONS_H

#include <vector>
#include <cassert>
#include <string>
#include <fstream>
#include "Island.h"
#include "STTtable.h"

double getLogGrowth(const Island& island) noexcept;

double extractSumOfRates(const Island& island);

Species findOldestSpecies(const std::vector<Species>& speciesVec);

STTtable mergeSTTtables(const std::vector<STTtable>& STTvec, const int& n_samples);

int howManyLineages(const std::vector<Species>& species);

std::vector<SpeciesID> whichMainAncestors(const std::vector<Species>& species);

void outputBranching(const Island fullIsland, const std::string& fileName);

void outputSTT(const STTtable& fullSTT, const std::string& fileName);

#endif //ARCHIDAISIE_FREE_FUNCTIONS_H
