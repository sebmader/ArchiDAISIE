//
// Created by Sebastian Mader on 05.09.2018.
//
// Review by Richel J.C. Bilderbeek
#ifndef ARCHIDAISIE_ISLAND_H
#define ARCHIDAISIE_ISLAND_H

#include <vector>
#include <cassert>
#include "DrawEvents.h"
#include "event_type.h"
#include "SpeciesID.h"
#include "Species.h"
#include "event_type.h"

//RJCB: indent comments correctly

class Island {
public:
    explicit Island(int k);

    int getCarryingCap() const noexcept;
    std::vector<SpeciesID> getSpeciesIDs() const noexcept;
    int getNSpecies() const noexcept;
    double returnLogGrowth();        // returns the logistic growth term (1-n/K)
        // RJCB: should be free function
        // (public interface can already supply all info for this calculation)
    std::vector<double> getLocalRates() { return mLocalRates; }
        // RJCB: AFAICS, this should not be cached, it is too early for an unproven
        // speed optimization yet. Calculate this when needed

    void addSpecies(const Species &);  // adds new species to species vector
    void deleteSpecies(const SpeciesID&);  // deletes species from species vector
    // TODO RJCB: why int? 'const Species&' is most natural, else 'const SpeciesID&' would be most natural
                                     // UNORDERD -> swap with last and pop_back
    int findPos(const SpeciesID&) const;  // find the position of //RJCB: implementation detail, make private
                    // species (input) in island vector
    const Species& findSpecies(SpeciesID) const;  // find the position //RJCB: why int? 'const Species&' is most natural, else 'const SpeciesID&' would be most natural
                    // of species (input) in island vector
                    // MUST find species -> otherwise will fail !!!

    void calculateIslRates(
            const std::vector<double>& islandParameters,
            const int& n_mainlandSpecies, // RJCB: Use 'const int' instead of 'const int&'.
            const int& n_islands,
            const double& sumLogGrowthWOthisIsl
    );              // initialise/calculate rates and store them in EventRates vector
                    // gam_i, gam_m, lamb_cl, lamb_al, mu_l
                    // per island -> doesn't include global rates !!!
    double extractSumOfRates() const noexcept;  // return the per-island rates vector
            // RJCB: should be a free function

    event_type sampleLocalEvent(std::mt19937_64);
                    // in case a local event is drawn, sample island, event and species
                    // it happens to

    void immigrate(const SpeciesID&, double);  // mainland species immigrates to that island
    int drawMigDestinationIsland(
        int originIsland, //RJCB: doesn't Island have a const int ID?
        std::vector<double>& LogGrowthTerms, //RJCB: AFAICS, can be calculated when needed
        const double& initialMigrationRate,
        std::mt19937_64
    );
                    // draw island species migrates to
    void migrate(Species, const double& );  // migrating onto this island
            // RJCB: Use 'const double' instead of 'const double&'
    void speciateClado(const SpeciesID&, double, SpeciesID&);
                    // island species cladogenetically speciates
    void speciateAna(const SpeciesID&, SpeciesID&);
                    // island species anagenetically speciates
    void goExtinct(const SpeciesID&);  // island species goes extinct

    const std::vector<Species>& getSpecies() const { return mSpecies; }
    void printIsland();  // prints island vector of species to the screen
            // RJCB: use operator<< instead

    void addIsland(const Island&);     // add island to THIS island
            // RJCB: cannot add islands to islands in reality
    // ### CAUTION ### : illogical! -> make this a nonmember function


private:
    std::vector<Species> mSpecies;  // phylogeny vector of species on island
    std::vector<double> mLocalRates;  // vector of rates
                    // for events PER ISLAND (5 rates: gamI, gamM, lambC, lambA, mu)
            // RJCB: AFAICS, this should not be cached, it is too early for an unproven
            // speed optimization yet. Calculate this when needed
    int mK; // Carrying capacity (should be const one day)
                    // for now: mIslandk = mK / n_islands
};

#endif // ARCHIDAISIE_ISLAND_H
