//
// Created by Sebastian Mader on 06.09.2018.
//

#ifndef ARCHIDAISIE_ARCHIPELAGO_H
#define ARCHIDAISIE_ARCHIPELAGO_H

#include <random>
#include <vector>
#include <iostream>
#include "DrawEvents.h"
#include "SpeciesID.h"
#include "event_type.h"
#include "Species.h"
#include "Island.h"

class Archipelago {         // class for the whole archipelago
public:
    explicit Archipelago(const int & = 0, const int & = 0);
        // constructor of archipelago based on number of islands and Island(!)-wide K
//    void updateAliveSpec();    // updating the archipelago-wide vector of extant species

    int getNIslands() const noexcept;
    int getNSpecies() const;
    int getNSpeciesID();
    std::vector<SpeciesID> getSpeciesIDs();
    std::vector<SpeciesID> getGlobalSpeciesIDs() const;
    const std::vector<Island> & getIslands() const { return mIslands; }
    int getCarryingCap() const noexcept { return mK; }
    std::vector<double> getGlobalRates() const noexcept { return mGlobalRates; }
    int getNColonisations() const noexcept { return mNColonisations; }

    bool isGlobal(const SpeciesID& speciesID) const;
    int whereIsSpecies(const Species& species) const;
    bool hasSpecies(const SpeciesID& speciesID) const;
    std::vector<int> findIsl(const SpeciesID&) const;    // find the island(s) where
                                    // species (input) is within archipelago
    std::vector<Species> findIslSpecies(const SpeciesID& speciesID) const;
    std::vector<Species> findMostRecentSistersPops(const Species& species) const;
    void addSpecies(const Species&, int);  // adds species to island
    void updateNColonisations();

    void calculateAllRates(const std::vector<double>&,
            int n_mainlandSpecies, int n_islands);
                    // calculate per-island rates and global rates
                    // and save them in LocalRates and GlobalRates vector, resp.
                    // Also, output of sum of both global (.first) and local (.second) rates
    event_type sampleNextEvent(std::mt19937_64&);   // draw next event;
                                    // output -> {event(0-7), species(ID)(,island(0-i))}
                                    // if global event -> vector.size() = 2, if local -> size = 3

    void speciateGlobalClado(const SpeciesID&, std::mt19937_64&, SpeciesID& maxSpeciesID);
                    // island species cladogenetically speciates over all islands
                    // (one population gets replaced by two new species populations)
                    // with random separation of archipelago into two groups/populations
    void speciateGlobalAna(const SpeciesID&, SpeciesID& maxSpeciesID);
                    // island species collectively (on all islands) diverges
                    // from mainland ancestor
    void goGlobalExtinct(const SpeciesID&);  // island species
                    // goes extinct on all islands it occures on

    void correctSisterTaxaLocal(const SpeciesID&, int island);
    void correctSisterTaxaGlobal(const SpeciesID& extinctSpID);


    void doGlobalEvent(const event_type& globalEvent,
            const SpeciesID& speciesID,
            std::mt19937_64& prng,
            SpeciesID& maxSpeciesID);

    void doLocalEvent(const event_type& localEvent,
            const SpeciesID& speciesID,
            std::mt19937_64& prng,
            const double& time,
            SpeciesID& maxSpeciesID,
            int island,
            const double& iniMigrationRate);

    void doNextEvent(const event_type&,
            const double&,
            std::mt19937_64&,
            const double&,
            SpeciesID& maxSpeciesID,
            const std::vector<SpeciesID>& mainSpeciesIDs); // switch-statement
                    // that calls event functions to update the island vector

    void addArchi(const Archipelago&);  // add an island vector to this one / consolidate them

    Island makeArchiTo1Island() const;   // aggregate all islands in
                    // archipelago as it would be one and return it
    void printArchi();  // print archipelago islands to screen


private:
    std::vector<double> mGlobalRates;  // vector of rates for global events
                                       // (3 rates: lambC, lambA, mu)
    std::vector<Island> mIslands;  // vector of island phylogenies of whole archipelago
    int mK;    // carrying capacity of archipelago
    int mNColonisations;  // cumulative number of independant colonisations
};

#endif //ARCHIDAISIE_ARCHIPELAGO_H
