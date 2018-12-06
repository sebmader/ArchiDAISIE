//
// Created by Sebastian Mader on 05.09.2018.
//

#include "Island.h"
#include "free_functions.h"

using namespace std;


// ------ Island member functions ------

Island::Island(const int k) : mK{k}
{
    assert(k >= 0);
    mNColonisations = 0;
}

int Island::getCarryingCap() const noexcept
{
    return mK;
}

int Island::getNColonisations() const noexcept
{
    return mNColonisations;
}

int Island::getNSpecies() const noexcept
{
    return static_cast<int>(mSpecies.size());
}

void Island::printIsland() const
{
    for (auto &i : mSpecies)
        i.printSpec();
    cout << '\n';
}

int Island::findPos(const SpeciesID& speciesID) const
{   // find the position of certain species (input) in species vector
        // if not present: output = -1
    assert(speciesID.getSpeciesID() >= 0);
    const int n_species = static_cast<int>(mSpecies.size());
    for (int i = n_species-1; i >= 0; --i)
        if (mSpecies[i].getSpecID() == speciesID)
            return i;
    return -1;
}

const Species& Island::findSpecies(const SpeciesID speciesID) const
{
    assert(speciesID.getSpeciesID() >= 0);
    const int index = findPos(speciesID);
    assert(index >= 0 && index < getNSpecies());
    return mSpecies[index];
}

Species& Island::findRefSpecies(const SpeciesID& speciesID)
{
    assert(speciesID.getSpeciesID() >= 0);
    const int index = findPos(speciesID);
    assert(index >= 0 && index < getNSpecies());
    return mSpecies[index];
}

void Island::calculateIslRates(
        const std::vector<double>& islandParameters,
        const int n_mainlandSpecies,
        const int n_islands,
        const double& sumLogGrowthWOThis)
{   // calculates the per-island rates of events, outputs them (for immediate use)
    // and saves them within the island-class; input -> initial parameters
    // order of island parameter vector: gam_i, gam_m, lamb_cl, lamb_al, mu_l

    // calculate the per-island rates:
    const int n_alive = getNSpecies();

    const double logGrowth = getLogGrowth(*this);
    assert(logGrowth >= 0.0);
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
    assert(immigrationRate >= 0.0); assert(migrationRate >= 0.0);
    assert(localCladoRate >= 0.0); assert(localAnaRate >= 0.0);
    assert(localExtinctRate >= 0.0);

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
void Island::immigrate(const SpeciesID& speciesID, const double& time)
{   // immigration from the mainland to THIS island

    assert(speciesID.getSpeciesID() >= 0);
    assert(time >= 0.0);
    Species newSpecies(time, speciesID, speciesID, 'I', false, time, time,
            vector<char>());
    if (hasSpecies(speciesID)) {  // if extant -> re-immigration
        // ("re-setting the clock" (= BirthT))
        const int pos = findPos(speciesID);
        assert(pos < getNSpecies() && pos >= 0);
        mSpecies[pos] = newSpecies;
    }
    else {
        if (getNSpecies() + 1 > mK)
            throw logic_error("Immigration would make number of species"
                              " exceed carrying capacity.\n");
        addSpecies(newSpecies);
        ++mNColonisations;
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
            migrationRates[k] = 0.0;
            continue;
        }
        migrationRates[k] = (initialMigrationRate * n_speciesAlive * LogGrowthTerms[k])
            / (n_islands*n_islands - n_islands);
        assert(migrationRates[k] >= 0.0);
    }
    const int destinationIsland = drawDisEvent(migrationRates, prng);
    assert(destinationIsland != originIsland);
    assert(destinationIsland >= 0 && destinationIsland < n_islands);

    return destinationIsland;
}

void Island::migrate(const Species& oldSpecies, const double& time)
{
    assert(oldSpecies.isValid());
    assert(time >= 0.0);
    vector<char> newCladostates = oldSpecies.getCladoStates();
    newCladostates.push_back('b');
    Species newSpecies = Species(time, oldSpecies.getParID(),
            oldSpecies.getSpecID(), oldSpecies.getStatus(), true,
            oldSpecies.getAncestralBT(), oldSpecies.getColonisationT(),
            newCladostates);
    // inherit ancestral birthT of oldSpecies IF it has already migrated before

    const SpeciesID speciesID = oldSpecies.getSpecID();
    if(hasSpecies(speciesID)) {
        if (findSpecies(speciesID).getBirth() < oldSpecies.getBirth()) {
            // if re-migration: re-set clock if resident is younger
            // otherwise circular re-migration would make species younger and younger
            const int pos = findPos(speciesID);
            assert(pos >= 0 && pos < getNSpecies());
            mSpecies[pos].setBirth(time);
        }
    }
    else {  // if first migration: add species to island
        if (getNSpecies()+1 > mK)
            throw logic_error("Migration would make number of species"
                              " exceed carrying capacity.\n");
        addSpecies(newSpecies);
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
    assert(speciesID.getSpeciesID() >= 0);
    assert(time >= 0.0);

    // find species
    const Species oldSpecies = findSpecies(speciesID);
    assert(oldSpecies.isValid());

    // 2 new species:
    vector<char> newCladoStates = oldSpecies.getCladoStates();
    newCladoStates.push_back('a');
    Species newSpecies1 = Species(oldSpecies.getBirth(), oldSpecies.getParID(),
            maxSpeciesID.createNewSpeciesID(), 'C', false,
            oldSpecies.getBirth(), oldSpecies.getColonisationT(),
            newCladoStates);
    newCladoStates.back() = 'b';
    Species newSpecies2 = Species(time, oldSpecies.getParID(),
            maxSpeciesID.createNewSpeciesID(), 'C', false,
            time, oldSpecies.getColonisationT(), newCladoStates);

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
    assert(speciesID.getSpeciesID() >= 0);

    // find species
    const Species oldSpecies = findSpecies(speciesID);
    assert(oldSpecies.isValid());
    // new species
    Species newSpecies = Species(oldSpecies.getBirth(), oldSpecies.getParID(),
            maxSpeciesID.createNewSpeciesID(), 'A', false,
            oldSpecies.getBirth(), oldSpecies.getColonisationT(),
            oldSpecies.getCladoStates());
    // parent goes extinct & daugther gets added to island
    goExtinct(speciesID);
    addSpecies(newSpecies);
}

void Island::goExtinct(const SpeciesID& speciesID)
{   // species goes extinct
    if(!hasSpecies(speciesID))
        throw logic_error("Species does not exist on island.\n");
    assert(speciesID.getSpeciesID() >= 0);
    // find position of species
    const int pos = findPos(speciesID);
    assert(pos >= 0 && pos < getNSpecies());
    // remove species
    mSpecies.erase(mSpecies.begin() + pos);
}

void Island::addIsland(const Island& islNew)
{   // adds another island to THIS (for aggregating archipelagos)

    const vector<Species>& island2 = islNew.getSpecies();

    // add species vector to THIS island
    if (!island2.empty()) {
        mSpecies.reserve(mSpecies.size() + island2.size());
        mSpecies.insert(mSpecies.end(), island2.begin(), island2.end());

        // delete duplicates
        for (int j = 0; j < getNSpecies(); ++j) {
            for (int k = j + 1; k < getNSpecies(); ++k) {
                assert(mSpecies[j].isValid() && mSpecies[k].isValid());
                if (mSpecies[j].getSpecID() ==
                        mSpecies[k].getSpecID()) {
                    if (mSpecies[j].isImmigrant() && !mSpecies[j].hasMigrated()) {
                        if (mSpecies[k].isImmigrant() && !mSpecies[k].hasMigrated()) {
                                // if both "pure" immigrants (= not migrated)
                                // take the most recent -> re-immigration
                            if (mSpecies[j].getBirth() <= mSpecies[k].getBirth()) {
                                assert(k >= 0);
                                deleteSpecies(k);
                                --k;
                            }
                            else {
                                assert(j >= 0);
                                deleteSpecies(j);
                                --j;
                                break;  // break out of inner loop to start the outer at same j
                            }
                        }
                        else {  // if j is pure immigrant but k is not
                                    // delete the non-immigrant (= migrant)
                            assert(mSpecies[k].hasMigrated());
                            assert(mSpecies[k].getStatus()=='I');
                            assert(k >= 0);
                            deleteSpecies(k);
                            --k;
                        }
                    }
                    else {  // if j is not immigrant OR has migrated!
                                // but k is -> k stays
                        if (mSpecies[k].isImmigrant() && !mSpecies[k].hasMigrated()) {
                            assert(j >= 0);
                            deleteSpecies(j);
                            --j;
                            break;  // break out of inner loop to start the outer at same j
                        }
                        else {  // both not pure immigrants -> older one stays
                            if (mSpecies[j].getBirth() >= mSpecies[k].getBirth()) {
                                assert(k >= 0);
                                deleteSpecies(k);
                                --k;
                            }
                            else {
                                assert(j >= 0);
                                deleteSpecies(j);
                                --j;
                                break;  // break out of inner loop to start the outer at same j
                            }
                        }
                    }
                }
            }
        }
        // sort by birth time
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
        mNColonisations += islNew.getNColonisations();
    }
}

void Island::addSpecies(const Species& newSpecies)
{   // adds new species to island -> both species and alive species vector
    assert(newSpecies.isValid());
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
    assert(pos >= 0 && pos < getNSpecies());
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
    assert(pos < getNSpecies());
    return pos >= 0;
}
