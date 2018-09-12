//
// Created by Sebastian Mader on 05.09.2018.
//

#ifndef ARCHIDAISIE_ISLAND_H
#define ARCHIDAISIE_ISLAND_H

#include <vector>
#include <cassert>
#include "DrawEvents.h"
#include "Species.h"
//#include "Archipelago.h"

//using namespace std;

class Island {        // class for ONE island within archipelago
public:
    explicit Island(const int k) : mIK{k} {assert(k >= 0);} // island constructor based on island-wide K

    int getCarryingCap() const noexcept { return mIK; }
    int getNSpeciesAlive() const noexcept {return static_cast<int>(mvIslSpecAlive.size());} // returns number of species alive
    int getNAllSpecies() const noexcept { return this->mvIsland.size(); }
    double returnLogGrowth() { return 1 - static_cast<double>(mvIslSpecAlive.size()) / mIK;}   // returns the logistic growth term (1-n/K)

    void addSpecies(const Species& newSpecies);    // adds new species to species vector
    int findPos(const int &ID) const;    // find the position of certain species (input) in IslandPhylo vector
    const Species& findSpecies(const int species_id) const;    // find the position of certain species (input) in IslandPhylo vector
    const Species& returnSpecies(const int pos) const;   // returns specific species from species vector

    int createNewID();       // returns new species ID and maxID += 1

    double calculateIslRates(const std::vector<double> &, const int &, const int &, const double &);
    // initialise/calculate rates and store them in EventRates vector
    // gam_i, gam_m, lamb_cl, lamb_al, mu_l
    // per island -> doesn't include global rates !!!
    double extractSumIslRate() const noexcept;      // return the per-island rates vector

    std::vector<int> sampleLocalEvent(std::mt19937_64, const int&);   // in case a local event is drawn, sample island, event and species
    // it happens to

    void immigrate(const int&, double);                   // mainland species immigrates to that island
    int drawMigrationIsland(const int originIsland, std::vector<double>&, const double&, std::mt19937_64);                     // island species migrates to other island
    void speciateClado(const int&, double);               // island species cladogenetically speciates
    void speciateAna(const int&, double);                 // island species anagenetically speciates
    void goExtinct(const int&, double);                   // island species goes extinct

    const std::vector<Species>& returnIsland() const { return mvIsland; }    // return island vector of species
    const std::vector<int>& returnIslSpecAlive() const { return mvIslSpecAlive; }  // return extant species vector
    void printIsland();                 // prints island vector of species to the screen

    void addIsland(const Island &);     // add island to THIS island
    // ### CAUTION ### : illogical! -> make this a nonmember function


private:
    std::vector<Species> mvIsland;           // phylogeny vector of species on island
    std::vector<int> mvIslSpecAlive;            // vector of species' IDs of alive species (.size() = n species on island)
    std::vector<double> mvLocalRates;        // vector of rates for events PER ISLAND (5 rates: gamI, gamM, lambC, lambA, mu)
    int mIK; // Carrying capacity (should be const one day)
    // for now: mIk = mAK / iNumIslands
};

#endif // ARCHIDAISIE_ISLAND_H
