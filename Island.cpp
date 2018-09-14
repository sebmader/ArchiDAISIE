//
// Created by Sebastian Mader on 05.09.2018.
//

#include "Island.h"

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
    const int n_alive = getNSpeciesAlive();
    int n_immiSpecies = 0;    // number of immigrant species
                              // -> only ones that can undergo anagenesis
    for (auto& g : mIsland)
        if (g.readSpID() <= 1000)
            ++n_immiSpecies;

    // immigration to specific island:
    const double immigrationRate = max(0.0, islandParameters[0] * n_mainlandSpecies
            * (1 - static_cast<double>(n_alive) / mIslandK) / n_islands);
    // migration from this island to all other islands:
    const double migrationRate = max(0.0, (islandParameters[1] * n_alive
            * sumLogGrowthWOthisIsl) / n_islands*n_islands - n_islands);
        // ### CAUTION ### : all other islands: '/iNumIsl',
                // one-way: '/iNumIsl^2 - iNumIsl' ???
    // local cladogenesis:
    const double localCladoRate = max(0.0, islandParameters[2] * n_alive
            * (1 - (static_cast<double>(n_alive) / mIslandK)));
    // local anagenesis:
    const double localAnaRate = max(0.0, islandParameters[3] * n_immiSpecies);
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
        vector<int> aliveSpecies;
        for(auto& species : mIsland) {
            if (species.isExtant())
                aliveSpecies.push_back(species.readSpID());
        }
        const int pos = drawUniEvent(0,
                static_cast<int>(aliveSpecies.size()), prng);
        speciesID = aliveSpecies[pos];
    }
    // return event and specID
    vector<int> vHappening = { event, speciesID };

    return vHappening;
}

// local updates:
void Island::immigrate(const int& speciesID, double time)
{   // immigration from the mainland to THIS island
    const int iPos = findPos(speciesID);  // check if species is
                                          // already present on island
    Species newSpecies(time, speciesID, speciesID);
    if (iPos >= 0) {  // if present
        assert(iPos < static_cast<int>(mIsland.size()));
        if (!mIsland[iPos].isExtant()) {   // if extinct
            mIsland.push_back(newSpecies);
        //    mvIslSpecAlive.push_back(speciesID);
        }
        else  // if extant -> re-immigration
              // ("re-setting the clock" (= BirthT))
            mIsland[iPos] = newSpecies;
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
        if (k == originIsland) {     // for the island it migrates from: LogGrowth = 0
            // -> this event doesn't happen
            migrationRates[k] = 0;
            continue;
        }
        migrationRates[k] = max(0.0, (initialMigrationRate * n_speciesAlive * LogGrowthTerms[k])
            / n_islands*n_islands - n_islands);
    }
    const int destinationIsland = drawDisEvent(migrationRates, prng);
    assert(destinationIsland != originIsland);

    return destinationIsland;
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
    // parent goes extinct
    goExtinct(speciesID, time);
    addSpecies(newSpecies1);
    addSpecies(newSpecies2);
}

void Island::speciateAna(const int& speciesID, double time, SpeciesID& maxSpeciesID)
{   // anagenetic speciation: only immigrants can undergo anagenesis !!!
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

void Island::addIsland(const Island &islNew)
{   // adds another island to THIS (for aggregating archipelagos)
    // extract data frames from island that's to be added
    const vector<Species>& vSpecNew = islNew.returnIsland();
    // const vector<int>& vSpecAliveNew = islNew.returnIslSpecAlive();

    // add species vector to THIS island
    if (!vSpecNew.empty()) {
        // intermediate vector of species
        vector<Species> vTempSpec = mIsland;
        cerr << mIsland.size() << '\n';
        cerr << vSpecNew.size() << '\n';
        vTempSpec.reserve(mIsland.size() + vSpecNew.size());
        vTempSpec.insert(vTempSpec.end(), vSpecNew.begin(), vSpecNew.end());
        // mIsland.reserve(mIsland.size() + vSpecNew.size());   // preallocate memory
        // mIsland.insert(mIsland.end(), vSpecNew.begin(), vSpecNew.end());
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

const Species& Island::returnSpecies(const int pos) const
{
  assert(pos >= 0);
  assert(pos < static_cast<int>(mIsland.size()));
  return mIsland[pos];
}   // returns specific species from species vector

void Island::addSpecies(const Species& newSpecies)
{   // adds new species to island -> both species and alive species vector
    mIsland.push_back(newSpecies);

//    const int speciesID = newSpecies.readSpID();
//    mvIslSpecAlive.push_back(speciesID);
}

std::vector<int> Island::getIDsSpeciesAlive() const
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
