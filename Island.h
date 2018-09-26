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
    explicit Island(const int k) : mIslandK{k} {assert(k >= 0);}
                    // island constructor based on island-wide K

    int getCarryingCap() const noexcept { return mIslandK; }
    int getNSpeciesAlive() const; // returns number of species alive
    std::vector<int> getIDsSpeciesAlive() const;
    int getNAllSpecies() const noexcept { return static_cast<int>(mIsland.size()); }
    double returnLogGrowth();  // returns the logistic growth term (1-n/K)

    void addSpecies(const Species &);  // adds new species to species vector
    int findPos(const int &) const;  // find the position of
                    // species (input) in island vector
//    int findPosAlive(const int &ID) const;  // find position of species in
                    // AliveSpecies vector
    const Species& findSpecies(int) const;  // find the position
                    // of species (input) in island vector
                    // MUST find species -> otherwise will fail !!!

    double calculateIslRates(const std::vector<double>&, const int&,
                                          const int&, const double&);
                    // initialise/calculate rates and store them in EventRates vector
                    // gam_i, gam_m, lamb_cl, lamb_al, mu_l
                    // per island -> doesn't include global rates !!!
    double extractSumOfRates() const noexcept;  // return the per-island rates vector

    std::vector<int> sampleLocalEvent(std::mt19937_64, const int&);
                    // in case a local event is drawn, sample island, event and species
                    // it happens to

    void immigrate(const int&, double);  // mainland species immigrates to that island
    int drawMigDestinationIsland(int, std::vector<double>&, const double&, std::mt19937_64);
                    // draw island species migrates to
    void migrate(const Species&);  // migrating onto this island
    void speciateClado(const int&, double, SpeciesID &);
                    // island species cladogenetically speciates
    void speciateAna(const int&, double, SpeciesID &);
                    // island species anagenetically speciates
    void goExtinct(const int&, double);  // island species goes extinct

    const std::vector<Species>& returnIsland() const { return mIsland; }
//    const std::vector<int>& returnIslSpecAlive() const { return mvIslSpecAlive; }
    void printIsland();  // prints island vector of species to the screen

    void consolidateIslands(const Island&);     // add island to THIS island
    // ### CAUTION ### : illogical! -> make this a nonmember function

private:
    std::vector<Species> mIsland;  // phylogeny vector of species on island
//    std::vector<int> mvIslSpecAlive;  // vector of species' IDs of alive species
    std::vector<double> mLocalRates;  // vector of rates
                    // for events PER ISLAND (5 rates: gamI, gamM, lambC, lambA, mu)
    int mIslandK; // Carrying capacity (should be const one day)
                    // for now: mIslandk = mArchiK / n_islands
};

#endif // ARCHIDAISIE_ISLAND_H
