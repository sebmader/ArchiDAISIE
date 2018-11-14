//
// Created by Sebastian Mader on 07.11.2018.
//

#ifndef ARCHIDAISIE_ARCHIDAISIE_H
#define ARCHIDAISIE_ARCHIDAISIE_H

#include <cstdlib>
#include <cmath>
#include <cassert>
#include <vector>
#include <fstream>
#include <chrono>
#include <random>
#include <exception>
#include "SpeciesID.h"
#include "Species.h"
#include "Island.h"
#include "Archipelago.h"
#include "free_functions.h"
#include "STTtable.h"

Archipelago ArchiDAISIE_core(const double& islandAge,
        const std::vector<SpeciesID>& mainSpeciesIDs,
        const std::vector<double>& initialParameters,
        int archiCarryingCap,
        int n_islands,
        std::mt19937_64& prng,
        SpeciesID& maxSpeciesID,
        STTtable& STT);

std::vector<Island> ArchiDAISIE(const double& islandAge,
        int n_mainlandSpecies,
        std::vector<double> initialParameters,
        int n_islands,
        int replicates = 1,
        int n_timeSlicesSTT = 25);

#endif //ARCHIDAISIE_ARCHIDAISIE_H
