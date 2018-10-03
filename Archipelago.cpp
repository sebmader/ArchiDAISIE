//
// Created by Sebastian Mader on 06.09.2018.
//


#include "Archipelago.h"

using namespace std;

Archipelago::Archipelago(const int &n_islands, const int &islCarryingCap)
    // constructor of archipelago based on number
{   // of islands and island-wide K
        assert(n_islands >= 0);
        assert(islCarryingCap >= 0);
        mIslands = vector<Island>((unsigned) n_islands,
                Island(islCarryingCap));
        mK = islCarryingCap * n_islands;
    try {
            if (n_islands == 0 || islCarryingCap == 0)
                throw string("You are creating an archipelago without islands"
                             " and/or with a carrying capacity of zero.\n");
        }
    catch(string &message) {
        clog << "Warning: " << message << '\n';
    }
}

vector<int> Archipelago::findIsl(const SpeciesID& speciesID) const
// find the island(s) where certain species (input) is within archipelago
// returns vector with island IDs (position in mIslands vector)
{
    vector<int> locations;
    const int n_islands = getNIslands();
    for (int i = 0; i < n_islands; ++i) {
        if (mIslands[i].hasSpecies(speciesID))
            locations.push_back(i);
    }
    return locations;
}

void Archipelago::calculateAllRates(
        const vector<double>& initialParameters,
        const int& n_mainlandSpecies,
        const int& n_islands)
{   // the rates for each event (local/global) are calculated and saved
    // order of output (vector.size()=2): first -> global, second -> local
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
        // for giving it to the island rates calculation function

    // count extant GLOBAL species on archipelago, for global rates
    vector<SpeciesID> globalSpeciesIDs = getGlobalSpeciesIDs();
    int n_globalSpecies = static_cast<int>(globalSpeciesIDs.size());

    // count alive immigrant species; only ones that can undergo anagenesis
    int n_globalImmigrants = 0;
    for (auto& speciesID : globalSpeciesIDs)
        if (speciesID <= static_cast<SpeciesID>(n_mainlandSpecies))
            ++n_globalImmigrants;

    // count all alive species -> for carrying capacity
    const int aliveSpecies = getNSpecies();
    // calculate global rates:
    // global cladogenesis: carrying capacity of all species, not only global species
    const double globalCladoRate = max(0.0, initialGlobalClado * n_globalSpecies *
            (1 - (static_cast<double>(aliveSpecies)/mK)));
    // global anagenesis:
    const double globalAnaRate = max(0.0, initialGlobalAna * n_globalImmigrants);
    // global extinction:
    const double globalExtinctRate = max(0.0, initialGlobalExtinct * n_globalSpecies);

    mGlobalRates = { globalCladoRate, globalAnaRate, globalExtinctRate };  // save global rates

    // calculate local rates:
    double sumLogGrowth = 0.0;  // sum of logistic growth term for migration
    for (auto &j : mIslands) {
        sumLogGrowth += getLogGrowth(j);
    }
    // calculate and save local rates per island:
    for (auto &i : mIslands) {
        double sumLogGrowthWOThis = sumLogGrowth - getLogGrowth(i); // sum
                    // of log-growth of all except THIS island (i)
        i.calculateIslRates(initialIslandPars,
                n_mainlandSpecies, n_islands, sumLogGrowthWOThis);
    }
}

event_type Archipelago::sampleNextEvent(mt19937_64 prng)
{   // which event will happen next

    // draw local event
    const int n_islands = static_cast<int>(mIslands.size());
    int n_localEvents = static_cast<int>(mIslands[0].getLocalRates().size());

    // sum of event rates over all islands
    vector<double> sumRatesPerEvent(n_localEvents,0);
    for (int i = 0; i < n_islands; ++i) {
        vector<double> islRates = mIslands[i].getLocalRates();
        assert(static_cast<int>(islRates.size()) == n_localEvents);
        for (int j = 0; j < n_localEvents; ++j)
            sumRatesPerEvent[j] += islRates[j];
    }

    // add global rates to vector of sums of local rates
    sumRatesPerEvent.reserve(sumRatesPerEvent.size() + mGlobalRates.size());
    sumRatesPerEvent.insert(sumRatesPerEvent.end(), mGlobalRates.begin(), mGlobalRates.end());
    assert(sumRatesPerEvent.size() == mGlobalRates.size() + n_localEvents);

    return static_cast<event_type>(drawDisEvent(sumRatesPerEvent, prng));
}

// global events:

void Archipelago::speciateGlobalClado(const SpeciesID& speciesID,
        mt19937_64 prng,
        double time,
        SpeciesID& maxSpeciesID)
    // species (input) globally cladogenetically speciates
{   // -> archipelago population splits into two new species
    vector<int> onWhichIslands = findIsl(speciesID);  // vector with islandIDs
        // (position in mIslands) where species is present
    if (onWhichIslands.size() < 2)
        throw logic_error("Drawn species is present on less than 2 islands. "
                          "Something's wrong.. (global cladogenesis)\n");

    // two daughter species
    const Species sp = mIslands[onWhichIslands[0]].findSpecies(speciesID);
    const double birthT = sp.getBirth();
    Species spNew1(birthT, speciesID, maxSpeciesID.createNewSpeciesID(), 'C');
    Species spNew2(time, speciesID, maxSpeciesID.createNewSpeciesID(), 'C');

    // draw where to split the archipelago:
        // 0 to i-2 -> split after the island number drawn
    int n_inhabitedIslands = static_cast<int>(onWhichIslands.size());
    assert(n_inhabitedIslands >= 2);
    const int split = onWhichIslands[drawUniEvent(0, n_inhabitedIslands-2, prng)];

    // update data frame
    for (auto& isl : onWhichIslands) {
        mIslands[isl].goExtinct(speciesID);
        if (isl <= split)   // split: after island with position equal to 'split'
            mIslands[isl].addSpecies(spNew1);
        else
            mIslands[isl].addSpecies(spNew2);
    }
}

void Archipelago::speciateGlobalAna(const SpeciesID& speciesID, SpeciesID& maxSpeciesID)
{   // species (input) globally anagenetically speciates
        // -> whole archipelago population diverges from mainland sp
        // can only happen to immigrant species
    vector<int> onWhichIslands = findIsl(speciesID);  // vector with islandIDs
                        // (position in mIslands) where species is present
    if (onWhichIslands.size() < 2)
        throw logic_error("Drawn species is present on less than 2 islands. "
                          "Something's wrong.. (global anagenesis)\n");
    // daughter species
    const Species sp = mIslands[onWhichIslands[0]].findSpecies(speciesID);
    const double birthT = sp.getBirth();
    Species spNew(birthT, speciesID, maxSpeciesID.createNewSpeciesID(), 'A');
    for (auto& iIsl : onWhichIslands) {
        mIslands[iIsl].goExtinct(speciesID);
        mIslands[iIsl].addSpecies(spNew);
    }
}

void Archipelago::goGlobalExtinct(const SpeciesID& speciesID)
{   // one species (input) goes exinct on all islands it inhabits
    vector<int> onWhichIslands = findIsl(speciesID);
    if (onWhichIslands.size() < 2)
        throw logic_error("Drawn species is present on less than 2 islands. "
                          "Something's wrong.. (global extinction)\n");
    for (auto& Isl : onWhichIslands) {
        mIslands[Isl].goExtinct(speciesID);
    }
}

void Archipelago::doGlobalEvent(const event_type globalEvent,
        const SpeciesID speciesID,
        mt19937_64 prng,
        const double& time,
        SpeciesID& maxSpeciesID)
{
    switch(globalEvent) {
    case event_type::global_cladogenesis:
        assert(getEventInt(globalEvent) == 5);
        speciateGlobalClado(speciesID, prng, time, maxSpeciesID);
        break;
    case event_type::global_anagenesis:
        assert(getEventInt(globalEvent) == 6);
        speciateGlobalAna(speciesID, maxSpeciesID);
        break;
    case event_type::global_extinction:
        assert(getEventInt(globalEvent) == 7);
        goGlobalExtinct(speciesID);
        break;
    default:
        throw logic_error("Event is not global.\n");
    }
}

void Archipelago::doLocalEvent(const event_type localEvent,
        const SpeciesID speciesID,
        std::mt19937_64 prng,
        const double& time,
        SpeciesID& maxSpeciesID,
        const int island,
        const double& iniMigrationRate)
{
    switch(localEvent)
    {
    case event_type::local_immigration:
    {
        assert(getEventInt(localEvent) == 0);
        mIslands[island].immigrate(speciesID, time);
        break;
    }
    case event_type::local_migration:
    {
        assert(getEventInt(localEvent) == 1);
        // save the logarithmic growth terms for all islands per island
        const int n_islands = getNIslands();
        vector<double> vLogs(n_islands, 0);
        for (int j = 0; j < n_islands; ++j) {
            vLogs[j] = getLogGrowth(mIslands[j]);
        }
        int destinationIsland = mIslands[island].drawMigDestinationIsland(
                island, vLogs, iniMigrationRate, prng);
        // output: position of island in mIslands where the species
        // migrates to
        Species newSpecies = mIslands[island].findSpecies(speciesID);
        mIslands[destinationIsland].migrate(newSpecies, time);
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
        mIslands[island].goExtinct(speciesID);
        break;
    default:
        throw logic_error("Event is not local.\n");
    }
}

void Archipelago::doNextEvent(const event_type nextEvent,
        const double& initialMigrationRate,
        mt19937_64 prng,
        double time,
        SpeciesID& maxSpeciesID)
{   // updates data frames; based on output of sampleNextEvent (event_type)
    // order of parameter indexes (Event):
        // immigration (0), migration (1), localClado(2), localAna (3),
        // localExtinct (4), globalClado(5), GlobalAna(6), GlobalExtinct(7)

    if (is_global(nextEvent)) {
        // sample global species
        SpeciesID speciesID = drawUniEvent(getGlobalSpeciesIDs(),prng);
        doGlobalEvent(nextEvent, speciesID, prng, time, maxSpeciesID);
    }
    else if (is_local(nextEvent)) {
        // sample island:
        const int n_islands = getNIslands();
        vector<double> eventRatePerIsland(n_islands, 0);
        for (int i = 0; i < n_islands; ++i) {
            vector<double> localRates = mIslands[i].getLocalRates();
            eventRatePerIsland[i] = localRates[getEventInt(nextEvent)];
        }
        const int isl = drawDisEvent(eventRatePerIsland, prng);
        // sample species:
        vector<SpeciesID> aliveSpecies = mIslands[isl].getSpeciesIDs();
        SpeciesID speciesID = drawUniEvent(aliveSpecies, prng);
        doLocalEvent(nextEvent, speciesID, prng, time, maxSpeciesID, isl, initialMigrationRate);
    }
    else
        throw logic_error("Next event is neither global nor local.. Something is wrong.\n");
}

void Archipelago::addArchi(const Archipelago &newArchi)
{   // adds one archipelago data frame (mIslands) to this one
    // important for putting together the 1-coloniser-archipelagos
        // -> means, there are no duplicates ?!

    const vector<Island> &addArch = newArchi.getIslands();
    if(addArch.size() != mIslands.size())
        throw logic_error("You cannot merge archipelagos of different size.\n");

    // consolidate single islands together
    const int n_islands = getNIslands();
    for (int i = 0; i < n_islands; ++i) {
        if (!addArch[i].getSpecies().empty()) { // ### CAUTION ### : does this make sense ??
            mIslands[i].addIsland(addArch[i]);
        }
    }
}

vector<Species> Archipelago::makeArchiTo1Island() const
{   // put islands of archipelago together in one Island vector,
        // delete duplicates, sort by time of birth/immigration

    vector<Species> archiToIsland;
    // add all island species vectors together
    for(auto& island : mIslands) {
        if (!island.getSpecies().empty()) {
            const vector<Species>& vTmp = island.getSpecies();
            archiToIsland.reserve(archiToIsland.size() + vTmp.size());
            archiToIsland.insert(archiToIsland.begin(),
                    vTmp.begin(), vTmp.end());
        }
    }

    // delete duplicates; ### CAUTION ### : what birth time ?!
    for (int j = 0; j < static_cast<int>(archiToIsland.size()); ++j) {
        for (int k = j + 1; k < static_cast<int>(archiToIsland.size()); ++k)
            if (archiToIsland[j].getSpecID() ==
                    archiToIsland[k].getSpecID()) { // TODO
                // ExtinctTime: take the extant one, or the later extinction
                // BirthTime: take the oldest birth time (initial colonisation) or
                    // the latest re-immigration time.. ### CAUTION ### : How??
                archiToIsland[k] = archiToIsland.back();
                archiToIsland.pop_back();
                --k;
            }
    }
    // sort by birth time
    const int newVecSize = static_cast<int>(archiToIsland.size());
    for (int l = 0; l < newVecSize - 1; ++l) {
        for (int m = l + 1; m < newVecSize; ++m) {
            if(archiToIsland[l].getBirth() <archiToIsland[m].getBirth()) {
                const Species tmpSp = archiToIsland[m];
                archiToIsland[m] = archiToIsland[l];
                archiToIsland[l] = tmpSp;
            }
        }
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
             << "SpID" << '\t' << "Stac" << '\n';
        z.printIsland();
    }
}

int Archipelago::getNSpecies()
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
        for (int k = j + 1; k < static_cast<int>(aliveSpecies.size()); ++k)
            if (aliveSpecies[j] ==
                    aliveSpecies[k]) {
                aliveSpecies[k] = aliveSpecies.back();
                aliveSpecies.pop_back();
                --k;
            }
    }
    return aliveSpecies;
}

vector<SpeciesID> Archipelago::getGlobalSpeciesIDs() const
{
    vector<SpeciesID> aliveSpecies;
    for(auto& island : mIslands) {
        vector<SpeciesID> tmpIDs = island.getSpeciesIDs();
        if (tmpIDs.empty())
            continue;
        aliveSpecies.reserve(aliveSpecies.size() + tmpIDs.size());
        aliveSpecies.insert(aliveSpecies.end(), tmpIDs.begin(), tmpIDs.end());
    }
    // save duplicate species (= global species):
    vector<SpeciesID> aliveGlobalSpecies;
    const int vecSize = static_cast<int>(aliveSpecies.size());
    for (int j = 0; j < vecSize - 1; ++j) {
        for (int k = j + 1; k < vecSize; ++k)
            if (aliveSpecies[j] ==
                    aliveSpecies[k]) {
                aliveGlobalSpecies.push_back(aliveSpecies[j]);
                ++j;
            }
    }
    // if on more than 2 islands -> creates duplicates
    // remove duplicates:
    for (int j = 0; j < static_cast<int>(aliveGlobalSpecies.size()) - 1; ++j) {
        for (int k = j + 1; k < static_cast<int>(aliveGlobalSpecies.size()); ++k) {
            cout << aliveGlobalSpecies.size() << '\n';
            assert(k > j);
            if (aliveGlobalSpecies[j] == aliveGlobalSpecies[k]) {
                aliveGlobalSpecies[k] = aliveGlobalSpecies.back();
                aliveGlobalSpecies.pop_back();
                --k;
            }
        }
    }
    return aliveGlobalSpecies;
}

int Archipelago::getNIslands() const noexcept
{
    return static_cast<int>(mIslands.size());
}
