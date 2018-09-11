//
// Created by Bastophiles on 06.09.2018.
//

#ifndef ARCHIDAISIE_ARCHIPELAGO_H
#define ARCHIDAISIE_ARCHIPELAGO_H

#include <random>
#include <vector>
#include <iostream>
#include <random>
#include "DrawEvents.h"
#include "Species.h"
#include "Island.h"

class Archipelago {         // class for the whole archipelago
public:
    Archipelago(const int & = 0, const int & = 0);  // constructor of archipelago based on number of islands and
    // archipelago-wide K
    void updateAliveSpec();             // updating the archipelago-wide vector of extant species

    std::vector<int> findIsl(const int &) const;    // find the island(s) where certain species (input) is within archipelago

    std::vector<double> calculateAllRates(const std::vector<double> &, const int &iM, const int &iNumIsl);    // calculate per-island rates and global rates
    // and save them in LocalRates and GlobalRates vector, resp.
    // Also, output of sum of both global (.first) and local (.second) rates
    std::vector<int> sampleNextEvent(const std::vector<double> &, std::mt19937_64, const int &);   // draw next event; output -> {event(0-7), species(ID)(,island(0-i))}
    // if global event -> vector.size() = 2, if local -> size = 3
    int createNewID();       // returns new species ID and maxID += 1

    void speciateGlobalClado(const int&, std::mt19937_64, double dTime);         // island species cladogenetically speciates over all islands
    // (one population gets replaced by two new species populations)
    // with random separation of archipelago into two groups/populations
    void speciateGlobalAna(const int&, double dTime);           // island species collectively (on all islands) diverges from mainland ancestor
    void goGlobalExtinct(const int&, double dTime);             // island species goes extinct on all islands it occures on

    void updateArchi(
      const std::vector<int>&,
      const double&,
      std::mt19937_64,
      double
    ); // switch-statement that calls event functions
    // updates the ArchiPhylo vector
    // LOCAL events indicated by 3 elements: event ([0]), species ([1]), island ([2])
    // GLOBAL events indicated by 2 elements: event ([0]), species ([1])

    void addArchi(const Archipelago&);  // add an ArchiPhylo to this one / consolidate them

    std::vector<Species> aggregateArchi() const;   // aggregate all islands in ArchiPylo vector as it would be one and return it
    const std::vector<Island> & returnArchi() const {return mvArchipel;}      // return OneIslaArchiPylo vector (and output to file ?)
    // ### CAUTION ###: only const& if assigning this output to a variable copies it
    void printArchi();

    static void setMaxID(const int &maxID) {Archipelago::mMaxSpID = maxID;}
    static int& returnMaxID() {return Archipelago::mMaxSpID;}
    static void incrementMaxID() {++Archipelago::mMaxSpID;}

private:
    std::vector<double> mvGlobalRates;           // vector of rates for global events (3 rates: lambC, lambA, mu)
    std::vector<Island> mvArchipel;            // vector of island phylogenies of whole archipelago
    std::vector<int> mvArchSpecAlive;        // vector of IDs of extant species
    // ### CAUTION ### : how to update? -> at the end of each event / within the event update functions
    // std::vector<Species> mvOneIslArchipel;    // aggregated vector of phylogeny of whole archipelago (like it'd be one island)
    static int mMaxSpID;    // the currently highest species ID
    // ### CAUTION ### : static! --> one shared member var for all archipelago objects initialised
    int mAK;    // carrying capacity of archipelago
};


#endif //ARCHIDAISIE_ARCHIPELAGO_H
