//
// Created by Sebastian Mader on 05.09.2018.
//
// Review by Richel J.C. Bilderbeek
#ifndef ARCHIDAISIE_ISLAND_H
#define ARCHIDAISIE_ISLAND_H

#include <vector>
#include <cassert>
#include "DrawEvents.h"
#include "SpeciesID.h"
#include "Species.h"

//RJCB: indent comments correctly

class Island {        // class for ONE island within archipelago //RJCB: redundant comment: an Island is an island is one island
public:
    explicit Island(const int k) : mIslandK{k} {assert(k >= 0);} //RJCB: move implementation to .cpp
                    // island constructor based on island-wide K //RJCB: indent comments correctly

    int getCarryingCap() const noexcept { return mIslandK; }
//    int getNSpecies() const; // returns number of species alive //RJCB: remove unused code
    std::vector<int> getSpeciesIDs() const; //RJCB: is noexcept, as I cannot imagine this not to work
    int getNSpecies() const noexcept { return static_cast<int>(mIsland.size()); }
    double returnLogGrowth();  // returns the logistic growth term (1-n/K) //RJCB: should be free function (public interface can already supply all info for this calculation)
    std::vector<double> getLocalRates() { return mLocalRates; } //RJCB: AFAICS, this should not be cached, it is too early for an unproven speed optimization yet. Calculate this when needed

    void addSpecies(const Species &);  // adds new species to species vector
    void deleteSpecies(const int&);  // deletes species from species vector //RJCB: why int? 'const Species&' is most natural, else 'const SpeciesID&' would be most natural
                                     // UNORDERD -> swap with last and pop_back
    int findPos(const int &) const;  // find the position of //RJCB: implementation detail, make private
                    // species (input) in island vector
//    int findPosAlive(const int &ID) const;  // find position of species in //RJCB: remove unused code
                    // AliveSpecies vector //RJCB: remove unused code
    const Species& findSpecies(int) const;  // find the position //RJCB: why int? 'const Species&' is most natural, else 'const SpeciesID&' would be most natural
                    // of species (input) in island vector
                    // MUST find species -> otherwise will fail !!!

    double calculateIslRates(
      const std::vector<double>& islandParameters, //RJCB: should this not be a const member of the island?
      const int& n_mainlandSpecies, // RJCB: Use 'const int' instead of 'const int&'.
      const int& n_islands,
      const double& sumLogGrowthWOthisIsl //RJCB: I think the vector of Species allows to calculate that when needed
    );
                    // initialise/calculate rates and store them in EventRates vector
                    // gam_i, gam_m, lamb_cl, lamb_al, mu_l
                    // per island -> doesn't include global rates !!!
    double extractSumOfRates() const noexcept;  // return the per-island rates vector //RJCB: should be a free function

    std::vector<int> sampleLocalEvent(std::mt19937_64, const int& n_mainlandSpecies); // RJCB: should return an 'event_type', namely the event that will happen. Or should be private, as this is an implementation detail
                    // in case a local event is drawn, sample island, event and species
                    // it happens to

    void immigrate(const int&, double);  // mainland species immigrates to that island //RJCB: why int? 'const Species&' is most natural, else 'const SpeciesID&' would be most natural
    int drawMigDestinationIsland(
        const int originIsland, //RJCB: doesn't Island have a const int ID?
        vector<double>& LogGrowthTerms, //RJCB: AFAICS, can be calculated when needed
        const double& initialMigrationRate,
      std::mt19937_64
    );
                    // draw island species migrates to
    void migrate(Species, const double& );  // migrating onto this island //RJCB: Use 'const double' instead of 'const double&'
    void speciateClado(const int&, double, SpeciesID &);
                    // island species cladogenetically speciates
    void speciateAna(const int&, SpeciesID&);
                    // island species anagenetically speciates
    void goExtinct(const int&);  // island species goes extinct

    const std::vector<Species>& returnIsland() const { return mIsland; } //RJCB: ???rename to 'GetSpecies', this naming makes no sense
//    const std::vector<int>& returnIslSpecAlive() const { return mvIslSpecAlive; } //RJCB: remove unused code
    void printIsland();  // prints island vector of species to the screen //RJCB: use operator<< instead

    void consolidateIslands(const Island&);     // add island to THIS island //RJCB: cannot add islands to islands in reality
    // ### CAUTION ### : illogical! -> make this a nonmember function

private:
    std::vector<Species> mIsland;  // phylogeny vector of species on island //RJCB: rename to 'mSpecies' as species are not islands
    std::vector<double> mLocalRates;  // vector of rates //RJCB: AFAICS, this should not be cached, it is too early for an unproven speed optimization yet. Calculate this when needed
                    // for events PER ISLAND (5 rates: gamI, gamM, lambC, lambA, mu)
    int mIslandK; // Carrying capacity (should be const one day) //RJCB: redundant naming: the K of an Island can be called mK
                    // for now: mIslandk = mArchiK / n_islands
};

#endif // ARCHIDAISIE_ISLAND_H
