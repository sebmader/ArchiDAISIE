//
// Created by Bastophiles on 06.09.2018.
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
        // constructor of archipelago based on number of islands and archipelago-wide K
//    void updateAliveSpec();    // updating the archipelago-wide vector of extant species

    int getNSpecies();
    std::vector<int> getSpeciesIDs();
    std::vector<int> getGlobalSpeciesIDs();
    std::vector<int> findIsl(const int &) const;    // find the island(s) where
                                    // species (input) is within archipelago

    std::vector<double> calculateAllRates(const std::vector<double> &,
                    const int &n_mainlandSpecies, const int &n_islands);
                    // calculate per-island rates and global rates
                    // and save them in LocalRates and GlobalRates vector, resp.
                    // Also, output of sum of both global (.first) and local (.second) rates
    std::vector<int> sampleNextEvent(const std::vector<double> &,
                                    std::mt19937_64, const int &);   // draw next event;
                                    // output -> {event(0-7), species(ID)(,island(0-i))}
                                    // if global event -> vector.size() = 2, if local -> size = 3

    void speciateGlobalClado(const int&, std::mt19937_64,
                    double time, SpeciesID& maxSpeciesID);  // island species
                    // cladogenetically speciates over all islands
                    // (one population gets replaced by two new species populations)
                    // with random separation of archipelago into two groups/populations
    void speciateGlobalAna(const int&, SpeciesID& maxSpeciesID);
                    // island species collectively (on all islands) diverges
                    // from mainland ancestor
    void goGlobalExtinct(const int&);  // island species
                    // goes extinct on all islands it occures on

    void doNextEvent(const std::vector<int>&, const double&,
            std::mt19937_64, double, SpeciesID& maxSpeciesID); // switch-statement
                    // that calls event functions updates the ArchiPhylo vector
                    // LOCAL events indicated by 3 elements: { event, species, island }
                    // GLOBAL events indicated by 2 elements: { event, species }

    void addArchi(const Archipelago&);  // add an ArchiPhylo to this one / consolidate them

    std::vector<Species> makeArchiTo1Island() const;   // aggregate all islands in
                    // archipelago as it would be one and return it
    const std::vector<Island> & returnArchi() const {return mArchipel;}
    void printArchi();  // print archipelago islands to screen

    int getCarryingCap() const noexcept { return mArchiK; }
    std::vector<double> getGlobalRates() const noexcept { return mGlobalRates; }

private:
    std::vector<double> mGlobalRates;  // vector of rates for global events
                                       // (3 rates: lambC, lambA, mu)
    std::vector<Island> mArchipel;  // vector of island phylogenies of whole archipelago
    int mArchiK;    // carrying capacity of archipelago
};

#endif //ARCHIDAISIE_ARCHIPELAGO_H
