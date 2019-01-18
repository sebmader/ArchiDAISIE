//
// Created by Sebastian Mader on 06.09.2018.
//


#include "Archipelago.h"
#include "free_functions.h"

using namespace std;

Archipelago::Archipelago(const int &n_islands, const int &islCarryingCap)
    // constructor of archipelago based on number
{   // of islands and island-wide K
        assert(n_islands >= 0);
        assert(islCarryingCap >= 0);
        mIslands = vector<Island>((unsigned) n_islands,
                Island(islCarryingCap));
        mK = islCarryingCap * n_islands;
        mNColonisations = 0;
    try {
            if (n_islands == 0 || islCarryingCap == 0)
                throw string("You are creating an archipelago without islands"
                             " and/or with a carrying capacity of zero.\n");
        }
    catch(string &message) {
        clog << "Warning: " << message << '\n';
    }
}

int Archipelago::getNSpeciesID()
{
    return static_cast<int>(getSpeciesIDs().size());
}

vector<SpeciesID> Archipelago::getSpeciesIDs()
{
    vector<SpeciesID> aliveSpecies;
    for(auto& island : mIslands) {
        vector<SpeciesID> tmpIDs = island.getSpeciesIDs();
        aliveSpecies.reserve(aliveSpecies.size() + tmpIDs.size());
        aliveSpecies.insert(aliveSpecies.end(), tmpIDs.begin(), tmpIDs.end());
    }
    // remove duplicates:
    for (int j = 0; j < static_cast<int>(aliveSpecies.size() - 1); ++j) {
        for (int k = j + 1; k < static_cast<int>(aliveSpecies.size()); ++k) {
            assert(k > j);
            assert(k >= 0 && k < static_cast<int>(aliveSpecies.size()));
            assert(j >= 0 && j < static_cast<int>(aliveSpecies.size()));
            if (aliveSpecies[j] == aliveSpecies[k]) {
                aliveSpecies[k] = aliveSpecies.back();
                aliveSpecies.pop_back();
                --k;
            }
        }
    }
    return aliveSpecies;
}

vector<SpeciesID> Archipelago::getGlobalSpeciesIDs() const
{
    vector<SpeciesID> aliveSpeciesIDs;
    for(auto& island : mIslands) {
        vector<SpeciesID> tmpIDs = island.getSpeciesIDs();
        if (tmpIDs.empty())
            continue;
        aliveSpeciesIDs.reserve(aliveSpeciesIDs.size() + tmpIDs.size());
        aliveSpeciesIDs.insert(aliveSpeciesIDs.end(), tmpIDs.begin(), tmpIDs.end());
    }
    // save duplicate species IDs (= global species):
    vector<SpeciesID> aliveGlobalSpeciesIDs;
    const int vecSize = static_cast<int>(aliveSpeciesIDs.size());
    for (int j = 0; j < vecSize - 1; ++j) {
        for (int k = j + 1; k < vecSize; ++k)
            if (aliveSpeciesIDs[j] ==
                    aliveSpeciesIDs[k]) {
                aliveGlobalSpeciesIDs.push_back(aliveSpeciesIDs[j]);
                break;
            }
    }
    // if on more than 2 islands -> creates duplicates
    // remove duplicates:
    for (int j = 0; j < static_cast<int>(aliveGlobalSpeciesIDs.size()) - 1; ++j) {
        for (int k = j + 1; k < static_cast<int>(aliveGlobalSpeciesIDs.size()); ++k) {
            assert(k > j);
            assert(k >= 0 && k < static_cast<int>(aliveGlobalSpeciesIDs.size()));
            assert(j >= 0 && j < static_cast<int>(aliveGlobalSpeciesIDs.size()));
            if (aliveGlobalSpeciesIDs[j] == aliveGlobalSpeciesIDs[k]) {
                aliveGlobalSpeciesIDs[k] = aliveGlobalSpeciesIDs.back();
                aliveGlobalSpeciesIDs.pop_back();
                --k;
            }
        }
    }
    return aliveGlobalSpeciesIDs;
}

int Archipelago::getNIslands() const noexcept
{
    return static_cast<int>(mIslands.size());
}

bool Archipelago::isGlobal(const SpeciesID& speciesID) const
{
    assert(speciesID.getSpeciesID() >= 0);
    return findIsl(speciesID).size() >= 2;
}

int Archipelago::whereIsSpecies(const Species& species) const
{
    assert(species.isValid());
    for (size_t i = 0; i < mIslands.size(); ++i) {
        vector<Species> islSpecies = mIslands[i].getSpecies();
        for (auto& sp : islSpecies) {
            if (sp == species)
                return static_cast<int>(i);
        }
    }
    assert(!"Should not get here! Species is not present on archipelago."); //!OCLINT
    return -1;
}

bool Archipelago::hasSpecies(const SpeciesID& speciesID) const
{
    assert(speciesID.getSpeciesID() >= 0);
    vector<int> onWhichIsls = findIsl(speciesID);
    return !onWhichIsls.empty();
}

int Archipelago::getNSpecies() const
{
    int n_species = 0;
    for (auto& isl : mIslands) {
        n_species += isl.getNSpecies();
    }
    return n_species;
}

vector<int> Archipelago::findIsl(const SpeciesID& speciesID) const
// find the island(s) where certain species (input) is within archipelago
// returns vector with island IDs (position in mIslands vector)
{
    assert(speciesID.getSpeciesID() >= 0);
    vector<int> locations;
    const int n_islands = getNIslands();
    for (int i = 0; i < n_islands; ++i) {
        if (mIslands[i].hasSpecies(speciesID))
            locations.push_back(i);
    }
    return locations;
}

vector<Species> Archipelago::findIslSpecies(const SpeciesID& speciesID) const
// find the island(s) where certain species (input) is within archipelago
// returns vector with island IDs (position in mIslands vector)
{
    assert(speciesID.getSpeciesID() >= 0);
    vector<Species> populations;
    for (auto& isl : mIslands) {
        if (isl.hasSpecies(speciesID))
            populations.push_back(isl.findSpecies(speciesID));
    }
    return populations;
}

vector<Species> Archipelago::findMostRecentSistersPops(const Species& species) const
{  // including populations of same species
    assert(species.isValid());
    vector<Species> sisters;
    for(auto& isl : mIslands) {
        const vector<Species>& islSpecies = isl.getSpecies();
        for(auto& sp : islSpecies) {
            if (species.isMostRecentSis(sp))
                sisters.push_back(sp);
        }
    }
    return sisters;
}

void Archipelago::addSpecies(const Species& species, const int island)
{
    assert(species.isValid());
    assert(island < getNIslands());
    mIslands[island].addSpecies(species);
}

void Archipelago::updateNColonisations()
{
    // update number of archi-colonisations
    int tmpNColo = 0;
    for (int i = 0; i < getNIslands(); ++i) {
        tmpNColo += mIslands[i].getNColonisations();
    }
    assert(tmpNColo >= mNColonisations);
    mNColonisations = tmpNColo;
}

void Archipelago::calculateAllRates(
        const vector<double>& initialParameters,
        const int n_mainlandSpecies,
        const int n_islands)
{   // the rates for each event (local/global) are calculated and saved in class object
    // order of parameters (input): gam_i (0), gam_m (1), lamb_cl (2),
        // lamb_al (3), mu_l (4), lamb_cg (5), lamb_ag (6), mu_g (7)

    assert(initialParameters.size() == 8);
    // assign initial parameter values to variables
    const double initialGlobalClado = initialParameters[5];
    const double initialGlobalAna = initialParameters[6];
    const double initialGlobalExtinct = initialParameters[7];
    const vector<double> initialIslandPars = { initialParameters[0],
               initialParameters[1], initialParameters[2],
               initialParameters[3], initialParameters[4] };
    // order of island parameter vector: gam_i, gam_m, lamb_cl, lamb_al, mu_l

    // count extant GLOBAL species on archipelago, for global rates
    vector<SpeciesID> globalSpeciesIDs = getGlobalSpeciesIDs();
    int n_globalSpecies = static_cast<int>(globalSpeciesIDs.size());

    // count global immigrant species; only ones that can undergo anagenesis
    int n_globalImmigrants = 0;
    for (auto& speciesID : globalSpeciesIDs)
        if (speciesID <= static_cast<SpeciesID>(n_mainlandSpecies))
            ++n_globalImmigrants;

    // count all alive species -> for carrying capacity
    const int aliveSpecies = getNSpeciesID();
    // calculate global rates:
    // global cladogenesis: carrying capacity of all species, not only global species
    const double globalCladoRate = max(0.0, initialGlobalClado * n_globalSpecies *
            (1.0 - (static_cast<double>(aliveSpecies)/mK)));
    // global anagenesis:
    const double globalAnaRate = max(0.0, initialGlobalAna * n_globalImmigrants);
    // global extinction:
    const double globalExtinctRate = max(0.0, initialGlobalExtinct * n_globalSpecies);

    mGlobalRates = { globalCladoRate, globalAnaRate, globalExtinctRate };  // save global rates

    // calculate local rates:
    double sumLogGrowth = 0.0;  // sum of logistic growth term for migration
    for (auto& j : mIslands) {
        sumLogGrowth += getLogGrowth(j);
    }
    // calculate and save local rates per island:
    for (auto &i : mIslands) {
        const double sumLogGrowthWOThis = sumLogGrowth - getLogGrowth(i); // sum
                    // of log-growth of all islands except THIS (i)
        i.calculateIslRates(initialIslandPars,
                n_mainlandSpecies, n_islands, sumLogGrowthWOThis);
    }
}

event_type Archipelago::sampleNextEvent(mt19937_64& prng)
{   // which event will happen next

    const int n_islands = getNIslands();
    const int n_localEvents = static_cast<int>(mIslands[0].getLocalRates().size());

    // sum of rates per local events over all islands
    vector<double> sumRatesPerEvent((unsigned)n_localEvents,0.0);
    for (int i = 0; i < n_islands; ++i) {
        const vector<double> islRates = mIslands[i].getLocalRates();
        assert(static_cast<int>(islRates.size()) == n_localEvents);
        for (int j = 0; j < n_localEvents; ++j)
            sumRatesPerEvent[j] += islRates[j];
    }

    // add global rates to vector of sums of local rates
    sumRatesPerEvent.reserve(sumRatesPerEvent.size() + mGlobalRates.size());
    sumRatesPerEvent.insert(sumRatesPerEvent.end(), mGlobalRates.begin(), mGlobalRates.end());
    assert(sumRatesPerEvent.size() == mGlobalRates.size() + n_localEvents);

    for (double i : sumRatesPerEvent) {
        cout << i << ", ";
    }
    cout << '\n';

    return static_cast<event_type>(drawDisEvent(sumRatesPerEvent, prng));
}

// global events:

void Archipelago::speciateGlobalClado(const SpeciesID& speciesID,
        mt19937_64& prng,
        SpeciesID& maxSpeciesID)
    // species (input) globally cladogenetically speciates
{   // -> archipelago population splits into two new species
    assert(speciesID.getSpeciesID() >= 0);
    assert(maxSpeciesID.getSpeciesID() >= 0);

    vector<int> onWhichIslands = findIsl(speciesID);  // vector with islandIDs
        // (position in mIslands) where species is present
    if (onWhichIslands.size() < 2)
        throw logic_error("Drawn species is present on less than 2 islands. "
                          "Something's wrong.. (global cladogenesis)\n");

    // draw where to split the archipelago:
        // 0 to i-2 -> split after the island number drawn
    int n_inhabitedIslands = static_cast<int>(onWhichIslands.size());
    const int split = onWhichIslands[drawUniEvent(0, n_inhabitedIslands-2, prng)];

    // daughter IDs
    SpeciesID newSpeciesID1 = maxSpeciesID.createNewSpeciesID();
    SpeciesID newSpeciesID2 = maxSpeciesID.createNewSpeciesID();
    // update data frame
    for (auto& isl : onWhichIslands) {
        // parent populations are different on different islands (esp. birth time)
        const Species oldSpecies = mIslands[isl].findSpecies(speciesID);
        assert(oldSpecies.isValid());
        if (isl <= split) {   // split: after island with position equal to 'split'
            Species newSpecies1(oldSpecies.getBirth(), oldSpecies.getParID(),
                    newSpeciesID1, 'C', false,
                    oldSpecies.getBirth(), oldSpecies.getColonisationT(),
                    oldSpecies.getCladoStates());
            mIslands[isl].addSpecies(newSpecies1);
        }
        else {
            Species newSpecies2(oldSpecies.getBirth(), oldSpecies.getParID(),
                    newSpeciesID2, 'C', false,
                    oldSpecies.getBirth(), oldSpecies.getColonisationT(),
                    oldSpecies.getCladoStates());
            mIslands[isl].addSpecies(newSpecies2);
        }
        mIslands[isl].goExtinct(speciesID);
    }
}

void Archipelago::speciateGlobalAna(const SpeciesID& speciesID, SpeciesID& maxSpeciesID)
{   // species (input) globally anagenetically speciates
        // -> whole archipelago population diverges from mainland sp
        // can only happen to immigrant species
    assert(speciesID.getSpeciesID() >= 0);
    assert(maxSpeciesID.getSpeciesID() >= 0);

    vector<int> onWhichIslands = findIsl(speciesID);  // vector with islandIDs
                        // (position in mIslands) where species is present
    if (onWhichIslands.size() < 2)
        throw logic_error("Drawn species is present on less than 2 islands. "
                          "Something's wrong.. (global anagenesis)\n");
    const SpeciesID newSpeciesID = maxSpeciesID.createNewSpeciesID();
    // daughter species
    for (auto& isl : onWhichIslands) {
        const Species oldSpecies = mIslands[isl].findSpecies(speciesID);
        assert(oldSpecies.isValid());
        assert(oldSpecies.isImmigrant());  // as it only effects immigrants
        Species newSpecies(oldSpecies.getBirth(), oldSpecies.getParID(),
                newSpeciesID, 'A', false,
                oldSpecies.getBirth(), oldSpecies.getColonisationT());
        mIslands[isl].goExtinct(speciesID);
        mIslands[isl].addSpecies(newSpecies);
    }
}

void Archipelago::goGlobalExtinct(const SpeciesID& speciesID)
{   // one species (input) goes exinct on all islands it inhabits
    assert(speciesID.getSpeciesID() >= 0);
    vector<int> onWhichIslands = findIsl(speciesID);
    if (onWhichIslands.size() < 2)
        throw logic_error("Drawn species is present on less than 2 islands. "
                          "Something's wrong.. (global extinction)\n");
    for (auto& Isl : onWhichIslands) {
        mIslands[Isl].goExtinct(speciesID);
    }
}

void Archipelago::doGlobalEvent(const event_type& globalEvent,
        const SpeciesID& speciesID,
        mt19937_64& prng,
        SpeciesID& maxSpeciesID)
{
    assert(speciesID.getSpeciesID() >= 0);
    assert(maxSpeciesID.getSpeciesID() >= 0);

    if (!isGlobal(speciesID))
        throw logic_error("Drawn species is present on less than 2 islands. "
                          "Something's wrong.. (doGlobalEvent)\n");
    switch(globalEvent) {
    case event_type::global_cladogenesis:
        assert(getEventInt(globalEvent) == 5);
        speciateGlobalClado(speciesID, prng, maxSpeciesID);
        break;
    case event_type::global_anagenesis:
        assert(getEventInt(globalEvent) == 6);
        speciateGlobalAna(speciesID, maxSpeciesID);
        break;
    case event_type::global_extinction:
        assert(getEventInt(globalEvent) == 7);
        correctSisterTaxaGlobal(speciesID);
        goGlobalExtinct(speciesID);
        break;
    default:
        throw logic_error("Event is not global.\n");
    }
}

void Archipelago::doLocalEvent(const event_type& localEvent,
        const SpeciesID& speciesID,
        mt19937_64& prng,
        const double& time,
        SpeciesID& maxSpeciesID,
        const int island,
        const double& iniMigrationRate)
{
    assert(time >= 0.0);
    assert(speciesID.getSpeciesID() >= 0);
    assert(maxSpeciesID.getSpeciesID() >= 0);
    assert(iniMigrationRate >= 0.0);
    assert(island >= 0);

    switch(localEvent)
    {
    case event_type::local_immigration:
    {
        assert(getEventInt(localEvent) == 0);
        if(mIslands[island].hasSpecies(speciesID)) {
            correctSisterTaxaLocal(speciesID, island);
        }
        mIslands[island].immigrate(speciesID, time);
        break;
    }
    case event_type::local_migration:
    {
        assert(getEventInt(localEvent) == 1);
        // save the logarithmic growth terms for all islands per island
        // to draw the island of destination
        const int n_islands = getNIslands();
        vector<double> vLogs((unsigned) n_islands, 0.0);
        for (int j = 0; j < n_islands; ++j) {
            vLogs[j] = getLogGrowth(mIslands[j]);
        }
        int destinationIsland = mIslands[island].drawMigDestinationIsland(
                island, vLogs, iniMigrationRate, prng);
            // output: position of island in mIslands where the species
            // migrates to

        Species oldSpecies = mIslands[island].findSpecies(speciesID);
        // if the same species but an older coloniser is present
        // replace it with the younger coloniser
        if(mIslands[destinationIsland].hasSpecies(speciesID)) {
            Species resident = mIslands[destinationIsland].findSpecies(speciesID);
            if(oldSpecies.getColonisationT() < resident.getColonisationT() ||
                    oldSpecies.getBirth() > resident.getBirth()) {
                // in case of re-migration: need to correct the former sister's
                // daughter states -> might change the whole phylogeny
                correctSisterTaxaLocal(speciesID,destinationIsland);
                mIslands[destinationIsland].goExtinct(speciesID);
            }
        }
        oldSpecies = mIslands[island].findSpecies(speciesID);
        // if species is not present on island of destination
        // old population gets the new clado state 'a'
        if(!mIslands[destinationIsland].hasSpecies(speciesID)) {
            Species& refOldSpecies = mIslands[island].findRefSpecies(speciesID);
            vector<char> newCladoStates = refOldSpecies.getCladoStates();
            newCladoStates.push_back('a');
            refOldSpecies.setCladoStates(newCladoStates);
            assert(mIslands[island].findSpecies(speciesID).getCladoStates().size()
                     == newCladoStates.size());
        }
        mIslands[destinationIsland].migrate(oldSpecies, time);
        break;
    }
    case event_type::local_cladogenesis:
        assert(getEventInt(localEvent) == 2);
        mIslands[island].speciateClado(speciesID, time, maxSpeciesID);
        break;
    case event_type::local_anagenesis:
        assert(getEventInt(localEvent) == 3);
        mIslands[island].speciateAna(speciesID, maxSpeciesID);
        break;
    case event_type::local_extinction:
        assert(getEventInt(localEvent) == 4);
        // correct sister species
        correctSisterTaxaLocal(speciesID,island);
        mIslands[island].goExtinct(speciesID);
        break;
    default:
        throw logic_error("Event is not local.\n");
    }
    updateNColonisations();
}

void Archipelago::doNextEvent(const event_type& nextEvent,
        const double& initialMigrationRate,
        mt19937_64& prng,
        const double& time,
        SpeciesID& maxSpeciesID,
        const vector<SpeciesID>& mainSpeciesIDs)
{   // updates data frames; based on output of sampleNextEvent (event_type)
    // order of parameter indexes (Event):
        // immigration (0), migration (1), localClado(2), localAna (3),
        // localExtinct (4), globalClado(5), GlobalAna(6), GlobalExtinct(7)
    assert(time >= 0.0);
    assert(maxSpeciesID.getSpeciesID() >= 0);
    assert(initialMigrationRate >= 0.0);
    assert(!mainSpeciesIDs.empty());
    if (is_global(nextEvent)) {
        // sample global species
        vector<SpeciesID> globalSpecies = getGlobalSpeciesIDs();
        if (globalSpecies.empty())
            throw logic_error("No global species exist on archipelago "
                              "but drawn event is global.\n");
        SpeciesID speciesID = drawUniSpeciesID(globalSpecies, prng);
        if (nextEvent == event_type::global_anagenesis) { // sample only from global immigrants
            vector<SpeciesID> globalImmigrants;
            for (auto& globalSpID : globalSpecies) {
                if (globalSpID <= mainSpeciesIDs.back()) {
                    // works with 1-immigrant as well because that's the only species that
                    // can be a non-endemic
                    globalImmigrants.push_back(globalSpID);
                }
            }
            speciesID = drawUniSpeciesID(globalImmigrants,prng);
        }
        doGlobalEvent(nextEvent, speciesID, prng, maxSpeciesID);
    }
    else if (is_local(nextEvent)) {
        // sample island based on per-island-rate of drawn event:
        vector<double> eventRatePerIsland(mIslands.size(), 0.0);
        for (size_t i = 0; i < mIslands.size(); ++i) {
            const vector<double> localRates = mIslands[i].getLocalRates();
            eventRatePerIsland[i] = localRates[getEventInt(nextEvent)];
        }
        const int isl = drawDisEvent(eventRatePerIsland, prng);
        // sample species:
        SpeciesID speciesID = SpeciesID(drawUniSpeciesID(mainSpeciesIDs, prng));
        if (getEventInt(nextEvent)) { // -> if not immigration
            const vector<SpeciesID> aliveSpecies = mIslands[isl].getSpeciesIDs();
            speciesID = drawUniSpeciesID(aliveSpecies, prng);
        }
        doLocalEvent(nextEvent, speciesID, prng, time,
                maxSpeciesID, isl, initialMigrationRate);
    }
}

void Archipelago::addArchi(const Archipelago &newArchi)
{   // adds one archipelago data frame (mIslands) to this one
    // important for putting together the 1-coloniser-archipelagos
        // -> TODO: means, there are no duplicates, right?!

    const vector<Island> &addArch = newArchi.getIslands();
    if(addArch.size() != mIslands.size())
        throw logic_error("You cannot merge archipelagos of different size.\n");

    // consolidate single islands together
    const int n_islands = getNIslands();
    for (int i = 0; i < n_islands; ++i) {
        if (!addArch[i].getSpecies().empty()) {
            mIslands[i].addIsland(addArch[i]);
        }
    }
    updateNColonisations();
}

Island Archipelago::makeArchiTo1Island() const
{   // put islands of archipelago together in one Island vector,
        // delete duplicates, sort by time of birth/immigration

    Island archiToIsland(mK);
    // add all island species vectors together
    for(auto& island : mIslands) {
        archiToIsland.addIsland(island);
    }

    return archiToIsland;
}

void Archipelago::printArchi()
{
    int i = 0;
    for (auto& z : mIslands) {
        ++i;
        cout << "Island " << i << '\n';
        cout << "BirT" << '\t' << "ParID" << '\t'
             << "SpID" << '\t' << "Status" << '\t'
             << "Mig" << '\t'  << "AncesT" << '\t'
             << "ColoT" << '\t' << "CladoStac" << '\n';
        z.printIsland();
    }
}

void Archipelago::correctSisterTaxaGlobal(const SpeciesID& extinctSpID)
{
    assert(extinctSpID.getSpeciesID() >= 0);
    assert(isGlobal(extinctSpID));
    // find oldest population of extinct species
    const vector<Species> extinctPops = findIslSpecies(extinctSpID);
    Species oldestExtinctPop = findOldestSpecies(extinctPops);
    const int oldestPopIsl = whereIsSpecies(oldestExtinctPop);
    // find real sister species (not other populations of extinct species):
    // while there are other populations of extinctSp in that vector
    // delete the last daughter state of them -> get to first branching of other species
    bool isSameSpec;
    do {
        const vector<Species> sistersSpecies = findMostRecentSistersPops(oldestExtinctPop);
        isSameSpec = false;
        for (auto& sis : sistersSpecies) {
            if (sis.getSpecID()==oldestExtinctPop.getSpecID()) {
                isSameSpec = true;
                vector<char> newDaughterStates = oldestExtinctPop.getCladoStates();
                newDaughterStates.pop_back();
                oldestExtinctPop.setCladoStates(newDaughterStates);
                break;
            }
        }
    } while (isSameSpec);

    // update extinct species in island object
    Species& refOldestPop = mIslands[oldestPopIsl].findRefSpecies(oldestExtinctPop.getSpecID());
    refOldestPop.setCladoStates(oldestExtinctPop.getCladoStates());

    // now the local correctSisterTaxa does what's needed TODO: right??
    correctSisterTaxaLocal(extinctSpID, whereIsSpecies(oldestExtinctPop));
}

void Archipelago::correctSisterTaxaLocal(const SpeciesID& extinctSpID, const int island)
{
    const Species extinctSp = mIslands[island].findSpecies(extinctSpID);
    vector<char> extDaughterStates = extinctSp.getCladoStates();
    if (!extDaughterStates.empty()) {
        // if extinct species underwent cladogenesis at some point
        // all descendant sisters loose respective daughter state
        const vector<Species> sistersPops = findMostRecentSistersPops(extinctSp);
        Species oldestSisPop = findOldestSpecies(sistersPops);
        const int posLastSpeciation = static_cast<int>(extinctSp.getCladoStates().size())-1;
        for (auto& sis : sistersPops) {
            const int isl = whereIsSpecies(sis);
            Species& refSis = mIslands[isl].findRefSpecies(sis.getSpecID());
            if (sis==oldestSisPop && extDaughterStates.back() == 'a') { // TODO: correct?
                  // if daughter 'a' goes extinct, respective sister/population 'b' inherits
                  // its birth time and has not migrated
                refSis.setBirth(extinctSp.getBirth());
                assert(mIslands[isl].findSpecies(oldestSisPop.getSpecID()).getBirth()
                        ==extinctSp.getBirth());
                refSis.setMigrated(false);
                assert(!mIslands[isl].findSpecies(oldestSisPop.getSpecID()).hasMigrated());
            }
            // and all most recent daughters loose resp. daughter state
            vector<char> newDaughterStates = refSis.getCladoStates();
            assert(posLastSpeciation < (int)newDaughterStates.size());
            newDaughterStates.erase(newDaughterStates.begin() + posLastSpeciation);
            refSis.setCladoStates(newDaughterStates);
            assert(mIslands[isl].findSpecies(sis.getSpecID()).getCladoStates().size()
                    == sis.getCladoStates().size()-1);
        }
    }
}