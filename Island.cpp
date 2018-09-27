//
// Created by Sebastian Mader on 05.09.2018.
//

#include "Island.h"
#include "Archipelago.h"

using namespace std;


// ------ Island member functions ------

void Island::printIsland()
{
    for (auto &i : mIsland)
        i.printSpec();
    cout << '\n';
}

int Island::findPos(const int &speciesID) const
{   // find the position of certain species (input) in species vector
        // if not present: output = -1
    const int n_species = static_cast<int>(mIsland.size());
    for (int i = n_species-1; i >= 0; --i)
        if (mIsland[i].readSpID() == speciesID)
            return i;
    return -1;
}

const Species& Island::findSpecies(const int speciesID) const
{
    const int index = findPos(speciesID);
    assert(index >= 0);
    return mIsland[index];
}

double Island::calculateIslRates(
        const std::vector<double>& islandParameters,
        const int& n_mainlandSpecies,
        const int& n_islands,
        const double& sumLogGrowthWOthisIsl)
{   // calculates the per-island rates of events, outputs them (for immidiate use)
    // and saves them within the island-class; input -> initial parameters
    // order of island parameter vector: gam_i, gam_m, lamb_cl, lamb_al, mu_l

    // calculate the per-island rates
    const int n_alive = getNSpecies();

    const double logGrowth = returnLogGrowth();
    // immigration to specific island:
    const double immigrationRate = max(0.0, islandParameters[0] * n_mainlandSpecies
            * logGrowth / n_islands);
    // migration from this island to all other islands:
    const double migrationRate = max(0.0, (islandParameters[1] * n_alive
            * sumLogGrowthWOthisIsl) / (n_islands*n_islands - n_islands));
        // ### CAUTION ### : all other islands: '/iNumIsl',
                // one-way: '/iNumIsl^2 - iNumIsl' ???
    // local cladogenesis:
    const double localCladoRate = max(0.0, islandParameters[2] * n_alive
            * logGrowth);
    // local anagenesis:
    const double localAnaRate = max(0.0, islandParameters[3] * n_alive);
    // local extinction:
    const double localExtinctRate = max(0.0, islandParameters[4] * n_alive);

    // save in rate vector
    mLocalRates = { immigrationRate, migrationRate, localCladoRate,
                                   localAnaRate, localExtinctRate };
    double sumPerIslRates = immigrationRate + migrationRate + localCladoRate
                                         + localAnaRate + localExtinctRate;
    return sumPerIslRates;
}

double Island::extractSumOfRates() const noexcept
{
    double sumRates = 0;
    for (double rate : mLocalRates)
        sumRates += rate;
    return sumRates;
}

vector<int> Island::sampleLocalEvent(mt19937_64 prng,
        const int &n_mainlandSpecies)
{   // samples local event and species it happens to on THIS island
    // draw event
    const int event = drawDisEvent(mLocalRates, prng);
    // draw species
    int speciesID = drawUniEvent(1, n_mainlandSpecies, prng); // in case of immigration
    if (event) {   // if not immigration (1-4) -> SpecID from extant island species
        vector<int> aliveSpecies = getSpeciesIDs();
        const int pos = drawUniEvent(0,
                static_cast<int>(aliveSpecies.size() - 1) , prng);
        speciesID = aliveSpecies[pos];
    }
    // return event and specID
    vector<int> happening = { event, speciesID };

    return happening;
}

// local updates:
void Island::immigrate(const int& speciesID, double time)
{   // immigration from the mainland to THIS island
    const int pos = findPos(speciesID);  // check if species is
                                          // already present on island
    Species newSpecies(time, speciesID, speciesID, 'I');
    if (pos >= 0) {  // if extant -> re-immigration
                     // ("re-setting the clock" (= BirthT))
        assert(pos < static_cast<int>(mIsland.size()));
        mIsland[pos] = newSpecies;
    }
    else {
        mIsland.push_back(newSpecies);
    }
}

int Island::drawMigDestinationIsland(
        const int originIsland,
        vector<double>& LogGrowthTerms,
        const double& initialMigrationRate,
        mt19937_64 prng)
{   // migration from THIS island to another; output: island of destination
    // draw island to which species migrates
        // -> initial migration rate as parameter !!
    const int n_islands = static_cast<int>(LogGrowthTerms.size());
    const int n_speciesAlive = getNSpecies();
    vector<double> migrationRates(n_islands);
    for (int k = 0; k < n_islands; ++k) {
        if (k == originIsland) {     // for the island it migrates from:
                        // migration rate = 0 -> this event doesn't happen
            migrationRates[k] = 0;
            continue;
        }
        migrationRates[k] = (initialMigrationRate * n_speciesAlive * LogGrowthTerms[k])
            / (n_islands*n_islands - n_islands);
    }
    const int destinationIsland = drawDisEvent(migrationRates, prng);
    assert(destinationIsland != originIsland);

    return destinationIsland;
}

void Island::migrate(Species newSpecies, const double& time)
{
    newSpecies.setBirth(time);  // set birth time of migrating species
    // to the time of migration
    newSpecies.setStatus('M');
    const int speciesID = newSpecies.readSpID();
    const int pos = findPos(speciesID);
    if(pos == -1)  // if first migration: add species to island
        addSpecies(newSpecies);
    else {  // else (if re-migration): re-set clock // TODO: correct?
        assert(pos >= 0 && pos < static_cast<int>(mIsland.size()));
        mIsland[pos] = newSpecies;
    }
}

void Island::speciateClado(const int& speciesID, double time,
        SpeciesID& maxSpeciesID)
{   // island species cladogenetically diverges

    // find species
    const Species oldSpecies = findSpecies(speciesID);

    // 2 new species:
    Species newSpecies1 = Species(oldSpecies.readBirth(), speciesID,
            maxSpeciesID.createNewSpeciesID(), 'C');
    Species newSpecies2 = Species(time, speciesID,
            maxSpeciesID.createNewSpeciesID(), 'C');

    // parent goes extinct and daughters are added
    goExtinct(speciesID);
    addSpecies(newSpecies1);
    addSpecies(newSpecies2);
}

void Island::speciateAna(const int& speciesID, SpeciesID& maxSpeciesID)
{   // local anagenetic speciation: DAISIE interprets this as cladogenesis (?!),
    // so all species can undergo it not only migrants -> equivalent to global (=DAISIE) anagenesis
    // find species
    const Species oldSpecies = findSpecies(speciesID);
    // new species
    const double birthT = oldSpecies.readBirth();
    Species newSpecies = Species(birthT, speciesID,
            maxSpeciesID.createNewSpeciesID(), 'A');
    // parent goes extinct & daugther gets added to island
    goExtinct(speciesID);
    addSpecies(newSpecies);
}

void Island::goExtinct(const int& speciesID)
{   // species goes extinct
    // find species
    const int pos = findPos(speciesID);
    assert(pos >= 0);
    // remove species
    mIsland.erase(mIsland.begin() + pos);
}

void Island::consolidateIslands(const Island& islNew)
{   // adds another island to THIS (for aggregating archipelagos)

    const vector<Species>& island2 = islNew.returnIsland();
    // const vector<int>& vSpecAliveNew = islNew.returnIslSpecAlive();

    // add species vector to THIS island
    if (!island2.empty()) {
        mIsland.reserve(mIsland.size() + island2.size());
        mIsland.insert(mIsland.end(), island2.begin(), island2.end());

        // delete duplicates; ### CAUTION ### : what birth time ?!
        for (int j = 0; j < static_cast<int>(mIsland.size()); ++j) {
            for (int k = j + 1; k < static_cast<int>(mIsland.size()); ++k) { ;
                if (mIsland[j].readSpID() ==
                        mIsland[k].readSpID()) { // TODO
                    if (mIsland[j].isImmigrant()) {
                        if (mIsland[k].isImmigrant()) {  // if both immigrants
                                // take the most recent -> re-immigration
                            if (mIsland[j].readBirth() <= mIsland[k].readBirth()) {
                                deleteSpecies(k);
                                --k;
                            }
                            else {
                                deleteSpecies(j);
                                --j;
                            }
                        }
                        else {  // if j is immigrant but k is not
                                    // delete the non-immigrant (= migrant?!)
                            assert(mIsland[k].isMigrant());
                            deleteSpecies(k);
                            --k;
                        }
                    }
                    else {  // if j is not immigrant
                        if (mIsland[k].isImmigrant()) {  // but k is -> k stays
                            assert(mIsland[j].isMigrant());
                            deleteSpecies(j);
                            --j;
                        }
                        else {  // both not immigrants -> older one stays
                            assert(mIsland[j].isMigrant() || mIsland[k].isMigrant());
                            if (mIsland[j].readBirth() >= mIsland[k].readBirth()) {
                                deleteSpecies(k);
                                --k;
                            }
                            else {
                                deleteSpecies(j);
                                --j;
                            }
                        }
                    }
                }
            }
        }
    }
}

void Island::addSpecies(const Species& newSpecies)
{   // adds new species to island -> both species and alive species vector
    mIsland.push_back(newSpecies);
}

vector<int> Island::getSpeciesIDs() const
{
    vector<int> aliveSpecies;
    for (auto& species : mIsland)
        aliveSpecies.push_back(species.readSpID());
    return aliveSpecies;
}

double Island::returnLogGrowth()
{
    double logGrowth = 1.0 - static_cast<double>(getNSpecies()) / mIslandK;
    return logGrowth;
}

void Island::deleteSpecies(const int& pos)
{
    mIsland[pos] = mIsland.back();
    mIsland.pop_back();
}
