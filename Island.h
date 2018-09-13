//
// Created by Sebastian Mader on 05.09.2018.
//

#ifndef ARCHIDAISIE_ISLAND_H
#define ARCHIDAISIE_ISLAND_H

#include <vector>
#include <cassert>
#include "DrawEvents.h"
#include "SpeciesID.h"
#include "Species.h"

class Island {        // class for ONE island within archipelago
public:
    explicit Island(const int k) : mIslandK{k} {assert(k >= 0);} // island constructor based on island-wide K

    int getCarryingCap() const noexcept { return mIslandK; }
    int getNSpeciesAlive() const; // returns number of species alive
    std::vector<int> getIDsSpeciesAlive() const;
    int getNAllSpecies() const noexcept { return this->mIsland.size(); }
    double returnLogGrowth() { return 1.0 - static_cast<double>(getNSpeciesAlive()) / mIslandK;}   // returns the logistic growth term (1-n/K)

    void addSpecies(const Species& newSpecies);    // adds new species to species vector
    int findPos(const int &speciesID) const;    // find the position of certain species (input) in IslandPhylo vector
//    int findPosAlive(const int &ID) const; // find position of species in AliveSpecies vector
    const Species& findSpecies(const int speciesID) const;    // find the position of certain species (input) in IslandPhylo vector
    const Species& returnSpecies(const int pos) const;   // returns specific species from species vector

    double calculateIslRates(const std::vector<double>&, const int&, const int&, const double&);
    // initialise/calculate rates and store them in EventRates vector
    // gam_i, gam_m, lamb_cl, lamb_al, mu_l
    // per island -> doesn't include global rates !!!
    double extractSumOfRates() const noexcept;      // return the per-island rates vector

    std::vector<int> sampleLocalEvent(std::mt19937_64, const int&);   // in case a local event is drawn, sample island, event and species
    // it happens to

    void immigrate(const int&, double);                   // mainland species immigrates to that island
    int drawMigDestinationIsland(const int, std::vector<double>&, const double&, std::mt19937_64);                     // island species migrates to other island
    void speciateClado(const int&, double, SpeciesID& maxSpeciesID);               // island species cladogenetically speciates
    void speciateAna(const int&, double, SpeciesID& maxSpeciesID);                 // island species anagenetically speciates
    void goExtinct(const int&, double);                   // island species goes extinct

    const std::vector<Species>& returnIsland() const { return mIsland; }    // return island vector of species
//    const std::vector<int>& returnIslSpecAlive() const { return mvIslSpecAlive; }  // return extant species vector
    void printIsland();                 // prints island vector of species to the screen

    void addIsland(const Island &);     // add island to THIS island
    // ### CAUTION ### : illogical! -> make this a nonmember function

private:
    std::vector<Species> mIsland;           // phylogeny vector of species on island
//    std::vector<int> mvIslSpecAlive;            // vector of species' IDs of alive species (.size() = n species on island)
    std::vector<double> mLocalRates;        // vector of rates for events PER ISLAND (5 rates: gamI, gamM, lambC, lambA, mu)
    int mIslandK; // Carrying capacity (should be const one day)
    // for now: mIk = mArchiK / iNumIslands
};

#endif // ARCHIDAISIE_ISLAND_H
