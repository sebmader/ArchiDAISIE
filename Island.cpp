//
// Created by Sebastian Mader on 05.09.2018.
//

#include "Island.h"
#include "Archipelago.h"

using namespace std;


// ------ Island member functions ------

Island::Island(const int k) : mK{k}
{
    assert(k >= 0);
}

int Island::getCarryingCap() const noexcept
{
    return mK;
}

int Island::getNSpecies() const noexcept
{
    return static_cast<int>(mSpecies.size());
}


void Island::printIsland()
{
    for (auto &i : mSpecies)
        i.printSpec();
    cout << '\n';
}

int Island::findPos(const SpeciesID& speciesID) const
{   // find the position of certain species (input) in species vector
        // if not present: output = -1
    const int n_species = static_cast<int>(mSpecies.size());
    for (int i = n_species-1; i >= 0; --i)
        if (mSpecies[i].getSpecID() == speciesID)
            return i;
    return -1;
}

const Species& Island::findSpecies(const SpeciesID speciesID) const
{
    const int index = findPos(speciesID);
    assert(index >= 0);
    return mSpecies[index];
}

void Island::calculateIslRates(
        const std::vector<double>& islandParameters,
        const int& n_mainlandSpecies,
        const int& n_islands,
        const double& sumLogGrowthWOThis)
{   // calculates the per-island rates of events, outputs them (for immediate use)
    // and saves them within the island-class; input -> initial parameters
    // order of island parameter vector: gam_i, gam_m, lamb_cl, lamb_al, mu_l

    // calculate the per-island rates:
    const int n_alive = getNSpecies();

    const double logGrowth = getLogGrowth(*this);
    // immigration to specific island:
    const double immigrationRate = max(0.0, islandParameters[0] * n_mainlandSpecies
            * logGrowth / n_islands);
    // migration from this island to all other islands:
    const double migrationRate = max(0.0, (islandParameters[1] * n_alive
            * sumLogGrowthWOThis) / (n_islands*n_islands - n_islands));
                // one-way ('/iNumIsl^2 - iNumIsl') times sum of log growths of all other islands
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
}

event_type Island::sampleLocalEvent(mt19937_64& prng)
{   // samples local event and species it happens to on THIS island
    // draw event
    return static_cast<event_type>(drawDisEvent(mLocalRates, prng));
}

// local updates:
void Island::immigrate(const SpeciesID& speciesID, double time)
{   // immigration from the mainland to THIS island

    Species newSpecies(time, speciesID, speciesID, 'I', time);
    if (hasSpecies(speciesID)) {  // if extant -> re-immigration
        // ("re-setting the clock" (= BirthT))
        const int pos = findPos(speciesID);
        assert(pos < getNSpecies());
        assert(pos >= 0);
        mSpecies[pos] = newSpecies;
    }
    else {
        if (getNSpecies() + 1 > mK)
            throw logic_error("Immigration would make number of species"
                              " exceed carrying capacity.\n");
        addSpecies(newSpecies);
    }
}

int Island::drawMigDestinationIsland(
        const int originIsland,
        vector<double>& LogGrowthTerms,
        const double& initialMigrationRate,
        mt19937_64& prng)
{   // migration from THIS island to another; output: island of destination
    // draw island to which species migrates
        // -> initial migration rate as parameter !!
    const int n_islands = static_cast<int>(LogGrowthTerms.size());
    const int n_speciesAlive = getNSpecies();
    vector<double> migrationRates((unsigned)n_islands);
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

void Island::migrate(const Species& oldSpecies, const double& time)
{
    Species newSpecies = Species(time, oldSpecies.getParID(),
            oldSpecies.getSpecID(), 'M', oldSpecies.getCladeBirthT(),
            oldSpecies.getCladoStates());
    // save birth time of migrant as clade birth IF it hasn't already migrated before
    if (oldSpecies.getBirth() == oldSpecies.getCladeBirthT())
        newSpecies.setCladeBirth(oldSpecies.getBirth());

    const SpeciesID speciesID = oldSpecies.getSpecID();
    if(!hasSpecies(speciesID)) {  // if first migration: add species to island
        if (getNSpecies() + 1 > mK)
            throw logic_error("Migration would make number of species"
                              " exceed carrying capacity.\n");
        addSpecies(newSpecies);
    }
    else {  // else (if re-migration): re-set clock // TODO: correct?
        const int pos = findPos(speciesID);
        assert(pos >= 0 && pos < getNSpecies());
        mSpecies[pos] = newSpecies;
    }
}

void Island::speciateClado(const SpeciesID& speciesID, const double& time,
        SpeciesID& maxSpeciesID)
{   // island species cladogenetically diverges
    if(!hasSpecies(speciesID))
        throw logic_error("Species does not exist on island.\n");
    if (getNSpecies() + 1 > mK)
        throw logic_error("Cladogenesis would make number of species"
                          " exceed carrying capacity.\n");
    // find species
    const Species oldSpecies = findSpecies(speciesID);

    // 2 new species:
    vector<char> newCladoStates = oldSpecies.getCladoStates();
    newCladoStates.push_back('a');
    Species newSpecies1 = Species(oldSpecies.getBirth(), oldSpecies.getParID(),
            maxSpeciesID.createNewSpeciesID(), 'C', oldSpecies.getCladeBirthT(),
            newCladoStates);
    newCladoStates.back() = 'b';
    Species newSpecies2 = Species(time, oldSpecies.getParID(),
            maxSpeciesID.createNewSpeciesID(), 'C', oldSpecies.getCladeBirthT(),
            newCladoStates);

    // parent goes extinct and daughters are added
    goExtinct(speciesID);
    addSpecies(newSpecies1);
    addSpecies(newSpecies2);
}

void Island::speciateAna(const SpeciesID& speciesID, SpeciesID& maxSpeciesID)
{   // local anagenetic speciation: DAISIE interprets this as cladogenesis (?!),
    // so all species can undergo it not only migrants
    // -> equivalent to global (=DAISIE) anagenesis
    if(!hasSpecies(speciesID))
        throw logic_error("Species does not exist on island.\n");

    // find species
    const Species oldSpecies = findSpecies(speciesID);
    // new species
    const double birthT = oldSpecies.getBirth();
    Species newSpecies = Species(birthT, oldSpecies.getParID(),
            maxSpeciesID.createNewSpeciesID(), 'A', oldSpecies.getCladeBirthT(),
            oldSpecies.getCladoStates());
    // parent goes extinct & daugther gets added to island
    goExtinct(speciesID);
    addSpecies(newSpecies);
}

void Island::goExtinct(const SpeciesID& speciesID)
{   // species goes extinct
    if(!hasSpecies(speciesID))
        throw logic_error("Species does not exist on island.\n");

    // find species
    const int pos = findPos(speciesID);
    assert(pos >= 0);
    // remove species
    mSpecies.erase(mSpecies.begin() + pos);
}

void Island::addIsland(const Island& islNew)
{   // adds another island to THIS (for aggregating archipelagos)

    const vector<Species>& island2 = islNew.getSpecies();
    // const vector<int>& vSpecAliveNew = islNew.returnIslSpecAlive();

    // add species vector to THIS island
    if (!island2.empty()) {
        mSpecies.reserve(mSpecies.size() + island2.size());
        mSpecies.insert(mSpecies.end(), island2.begin(), island2.end());

        // delete duplicates; ### CAUTION ### : what birth time ?!
        for (int j = 0; j < getNSpecies(); ++j) {
            for (int k = j + 1; k < getNSpecies(); ++k) { ;
                if (mSpecies[j].getSpecID() ==
                        mSpecies[k].getSpecID()) {
                    if (mSpecies[j].isImmigrant()) {
                        if (mSpecies[k].isImmigrant()) {  // if both immigrants
                                // take the most recent -> re-immigration
                            if (mSpecies[j].getBirth() <= mSpecies[k].getBirth()) {
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
                            assert(mSpecies[k].isMigrant());
                            deleteSpecies(k);
                            --k;
                        }
                    }
                    else {  // if j is not immigrant
                        if (mSpecies[k].isImmigrant()) {  // but k is -> k stays
                            assert(mSpecies[j].isMigrant());
                            deleteSpecies(j);
                            --j;
                        }
                        else {  // both not immigrants -> older one stays
                            assert(mSpecies[j].isMigrant() || mSpecies[k].isMigrant());
                            if (mSpecies[j].getBirth() >=mSpecies[k].getBirth()) {
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
        // sort birth time
        const int newVecSize = static_cast<int>(mSpecies.size());
        for (int l = 0; l < newVecSize - 1; ++l) {
            for (int m = l + 1; m < newVecSize; ++m) {
                if(mSpecies[l].getBirth() <mSpecies[m].getBirth()) {
                    const Species tmpSp = mSpecies[m];
                    mSpecies[m] = mSpecies[l];
                    mSpecies[l] = tmpSp;
                }
            }
        }
    }
}

void Island::addSpecies(const Species& newSpecies)
{   // adds new species to island -> both species and alive species vector
    mSpecies.push_back(newSpecies);
}

vector<SpeciesID> Island::getSpeciesIDs() const
{
    vector<SpeciesID> aliveSpecies;
    for (auto& species : mSpecies)
        aliveSpecies.push_back(species.getSpecID());
    return aliveSpecies;
}

void Island::deleteSpecies(const int& pos)
{
    mSpecies[pos] = mSpecies.back();
    mSpecies.pop_back();
}

std::vector<double> Island::getLocalRates() const noexcept
{
    return mLocalRates;
}

bool Island::hasSpecies(const SpeciesID& speciesID) const
{
    const int pos = findPos(speciesID);
    return pos >= 0;
}
