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
#include "STT.h"
#include "STTtable.h"
#include "free_functions.h"

Archipelago ArchiDAISIE_core(const double& islandAge,
        const int n_mainlandSpecies,
        const std::vector<double>& initialParameters,
        const int archiCarryingCap,
        const int n_islands,
        std::mt19937_64& prng,
        SpeciesID& maxSpeciesID,
        STTtable& stt);

std::vector<Island> ArchiDAISIE(const double& islandAge,
        int n_mainlandSpecies,
        std::vector<double> initialParameters,
        int n_islands,
        int replicates);

#endif //ARCHIDAISIE_ARCHIDAISIE_H
