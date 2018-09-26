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
/*
int Island::findPosAlive(const int &speciesID) const
{
    const int aliveSpecies = getNSpeciesAlive();
    for (int i = 0; i < aliveSpecies; ++i)
        if (mvIslSpecAlive[i] == speciesID)
            return i;
    return -1;
}
*/

const Species& Island::findSpecies(const int speciesID) const
{
    const int index = findPos(speciesID);
    assert(index >= 0);
    return mIsland[index];
}

int Island::getNSpeciesAlive() const
{
    int extantSpeciesCounter = 0;
    for (auto& species : mIsland) {
        if (species.isExtant())
            ++extantSpeciesCounter;
    }
    return extantSpeciesCounter;
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
    const vector<int> aliveSpecies = getIDsSpeciesAlive();
    const int n_alive = static_cast<int>(aliveSpecies.size());

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
        vector<int> aliveSpecies = getIDsSpeciesAlive();
        const int pos = drawUniEvent(0,
                static_cast<int>(aliveSpecies.size()) - 1, prng);
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
    Species newSpecies(time, speciesID, speciesID);
    if (pos >= 0) {  // if present
        assert(pos < static_cast<int>(mIsland.size()));
        if (!mIsland[pos].isExtant()) {   // if extinct
            mIsland.push_back(newSpecies);
        //    mvIslSpecAlive.push_back(speciesID);
        }
        else  // if extant -> re-immigration
              // ("re-setting the clock" (= BirthT))
            mIsland[pos] = newSpecies;
    }
    else {
        mIsland.push_back(newSpecies);
    //    mvIslSpecAlive.push_back(speciesID);
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
    const int n_speciesAlive = getNSpeciesAlive();
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

void Island::migrate(const Species& newSpecies)
{
    const int speciesID = newSpecies.readSpID();
    const int pos = findPos(speciesID);
    if(pos == -1)  // if first migration: add species to island
        addSpecies(newSpecies);
            // else (if re-migration): nothing happens
}

void Island::speciateClado(const int& speciesID, double time,
        SpeciesID& maxSpeciesID)
{   // island species cladogenetically diverges

    // find species
    const Species oldSpecies = findSpecies(speciesID);

    // 2 new species:
    Species newSpecies1 = Species(oldSpecies.readBirth(), speciesID,
            maxSpeciesID.createNewSpeciesID());
    Species newSpecies2 = Species(time, speciesID,
            maxSpeciesID.createNewSpeciesID());

    // parent goes extinct and daughters are added
    goExtinct(speciesID, time);
    addSpecies(newSpecies1);
    addSpecies(newSpecies2);
}

void Island::speciateAna(const int& speciesID, double time, SpeciesID& maxSpeciesID)
{   // local anagenetic speciation: DAISIE interprets this as cladogenesis (?!),
    // so all species can undergo it not only migrants -> equivalent to global (=DAISIE) anagenesis
    // find species
    const Species oldSpecies = findSpecies(speciesID);
    // new species
    const double birthT = oldSpecies.readBirth();
    Species newSpecies = Species(birthT, speciesID,
            maxSpeciesID.createNewSpeciesID());
    // parent goes extinct & daugther gets added to island
    goExtinct(speciesID, time);
    addSpecies(newSpecies);
}

void Island::goExtinct(const int& speciesID, double time)
{   // species goes extinct
    // find species
    const int pos = findPos(speciesID);
    assert(pos >= 0);
    mIsland[pos].goExtinct(time);
//    const int posAlive = findPosAlive(speciesID);
//    mvIslSpecAlive[posAlive] = mvIslSpecAlive.back();
//    mvIslSpecAlive.pop_back();
}

void Island::consolidateIslands(const Island& islNew)
{   // adds another island to THIS (for aggregating archipelagos)
    // extract data frames from island that's to be added

    const vector<Species>& vSpecNew = islNew.returnIsland();
    // const vector<int>& vSpecAliveNew = islNew.returnIslSpecAlive();

    // add species vector to THIS island
    if (!vSpecNew.empty()) {
        mIsland.reserve(mIsland.size() + vSpecNew.size());
        mIsland.insert(mIsland.end(), vSpecNew.begin(), vSpecNew.end());
        // mIsland.reserve(mIsland.size() + vSpecNew.size());   // preallocate memory
        // mIsland.insert(mIsland.end(), vSpecNew.begin(), vSpecNew.end());
    }

    // delete duplicates; ### CAUTION ### : what birth time ?!
    const int islandSize = static_cast<int>(mIsland.size());
    for (int j = 0; j < islandSize; ++j) {
        for (int k = j + 1; k < islandSize; ++k)
            if (mIsland[j].readSpID() ==
                    mIsland[k].readSpID()) { // TODO
                // ExtinctTime: take the extant one, or the later extinction
                // BirthTime: take the oldest birth time (initial colonisation) or
                // the latest re-immigration time.. ### CAUTION ### : How??
                mIsland[k] = mIsland.back();
                mIsland.pop_back();
                --k;
            }
    }
    /*
    // add alive species to THIS island
    if (!vSpecAliveNew.empty()) {
        vector<int> vTempAlive = mvIslSpecAlive;
        vTempAlive.reserve(mvIslSpecAlive.size() + vSpecAliveNew.size());
        vTempAlive.insert(vTempAlive.end(), vSpecAliveNew.begin(), vSpecAliveNew.end());
        // mvIslSpecAlive.reserve(mvIslSpecAlive.size() + vSpecAliveNew.size());
        // mvIslSpecAlive.insert(mvIslSpecAlive.end(), vSpecAliveNew.begin(), vSpecAliveNew.end());
    }
    */
}

void Island::addSpecies(const Species& newSpecies)
{   // adds new species to island -> both species and alive species vector
    mIsland.push_back(newSpecies);

//    const int speciesID = newSpecies.readSpID();
//    mvIslSpecAlive.push_back(speciesID);
}

vector<int> Island::getIDsSpeciesAlive() const
{
    vector<int> aliveSpecies;
    for (auto& species : mIsland)
        if (species.isExtant())
            aliveSpecies.push_back(species.readSpID());
    return aliveSpecies;
}

double Island::returnLogGrowth()
{
    double logGrowth = 1.0 - static_cast<double>(getNSpeciesAlive()) / mIslandK;
    return logGrowth;
}