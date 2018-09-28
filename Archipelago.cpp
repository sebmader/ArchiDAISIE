//
// Created by Sebastian Mader on 06.09.2018.
//


#include "Archipelago.h"

using namespace std;

Archipelago::Archipelago(const int &n_islands, const int &archiCarryingCap)
    // constructor of archipelago based on number
{   // of islands and archipelago-wide K
    assert(n_islands >= 0); assert(archiCarryingCap >= 0);
    if(archiCarryingCap % n_islands != 0)
        throw logic_error("The carrying capacity needs to be dividable"
                          " by the number of islands without remainder.\n");
    if (n_islands == 0 || archiCarryingCap == 0)
        throw string("You are creating an archipelago without islands"
                     " and/or with a carrying capacity of zero.\n");
    mArchipel = vector<Island>((unsigned) n_islands,
            Island(archiCarryingCap / n_islands));
    mArchiK = archiCarryingCap;
}

vector<int> Archipelago::findIsl(const SpeciesID& speciesID) const
// find the island(s) where certain species (input) is within archipelago
// returns vector with island IDs (spots in mArchipel vector)
// ### CAUTION ### : maybe vector<pairs<int, int> > with island ID and position?
{
    vector<int> locations;
    const int n_islands = static_cast<int>(mArchipel.size());
    for (int i = 0; i < n_islands; ++i) {
        if (mArchipel[i].findPos(speciesID) >= 0)
            locations.push_back(i);
    }
    return locations;
}

vector<double> Archipelago::calculateAllRates(
        const vector<double> &initialParameters,
        const int &n_mainlandSpecies,
        const int &n_islands)
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
            (1 - (static_cast<double>(aliveSpecies)/mArchiK)));
    // global anagenesis:
    const double globalAnaRate = max(0.0, initialGlobalAna * n_globalImmigrants);
    // global extinction:
    const double globalExtinctRate = max(0.0, initialGlobalExtinct * n_globalSpecies);
    // sum of global rates:
    const double sumGlobal = globalCladoRate + globalAnaRate + globalExtinctRate;
    mGlobalRates = { globalCladoRate, globalAnaRate, globalExtinctRate };

    // calculate local rates:
    // logistic growth term for migration
    double sumLogGrowth = 0.0;
    for (auto &j : mArchipel) {
        sumLogGrowth += j.returnLogGrowth(); // sums the
            // logistic growth terms of all islands together
    }
    // sum of local rates:
    double sumLocal = 0;
    for (auto &i : mArchipel) {
        double sumLogGrowthWOthisIsl = sumLogGrowth - i.returnLogGrowth(); // sum
                    // of log-growth of all except THIS island (i)
        sumLocal += i.calculateIslRates(initialIslandPars,
                n_mainlandSpecies, n_islands, sumLogGrowthWOthisIsl);
    }

    vector<double> sumsGloLoc = { sumGlobal, sumLocal };

    return sumsGloLoc;
}

event_type Archipelago::sampleNextEvent(const vector<double>& sumsGloLoc,
        mt19937_64 prng,
        const int& n_mainlandSpecies)
{   // which event will happen next

    // draw local event summed up over all islands
    const int n_islands = static_cast<int>(mArchipel.size());
    int n_localEvents = static_cast<int>(mArchipel[0].getLocalRates().size());

    vector<double> sumRatesPerEvent(n_localEvents,0);
    for (int i = 0; i < n_islands; ++i) {
        vector<double> islRates = mArchipel[i].getLocalRates();
        assert(static_cast<int>(islRates.size()) == n_localEvents);
        for (int j = 0; j < n_localEvents; ++j)
            sumRatesPerEvent[j] += islRates[j];
    }

    sumRatesPerEvent.reserve(sumRatesPerEvent.size() + mGlobalRates.size());
    sumRatesPerEvent.insert(sumRatesPerEvent.end(), mGlobalRates.begin(), mGlobalRates.end());
    auto nextEvent = static_cast<event_type>(drawDisEvent(sumRatesPerEvent, prng));

    return nextEvent;
}

// global events:

void Archipelago::speciateGlobalClado(const SpeciesID& speciesID,
        mt19937_64 prng,
        double time,
        SpeciesID& maxSpeciesID)
    // species (input) globally cladogenetically speciates
{   // -> archipelago population splits into two new species
    vector<int> onWhichIslands = findIsl(speciesID);  // vector with islandIDs
        // (position in mArchipel) where species is present
    if (onWhichIslands.size() < 2)
        assert(!"Drawn species is present on less than 2 islands. "
                          "Something's wrong.. (global cladogenesis)\n");

    // two daughter species
    const Species sp = mArchipel[onWhichIslands[0]].findSpecies(speciesID);
    const double birthT = sp.readBirth();
    Species spNew1(birthT, speciesID, maxSpeciesID.createNewSpeciesID(), 'C');
    Species spNew2(time, speciesID, maxSpeciesID.createNewSpeciesID(), 'C');

    // draw where to split the archipelago:
        // 0 to i-2 -> split after the island number drawn
    int n_inhabitedIslands = static_cast<int>(onWhichIslands.size());
    assert(n_inhabitedIslands >= 2);
    const int split = drawUniEvent(0, n_inhabitedIslands-2, prng);

    // update data frame
    for (auto& isl : onWhichIslands) {
        mArchipel[isl].goExtinct(speciesID);
        if (isl <= split)   // split: after island with position equal to 'split'
            mArchipel[isl].addSpecies(spNew1);
        else
            mArchipel[isl].addSpecies(spNew2);
    }
}

void Archipelago::speciateGlobalAna(const SpeciesID& speciesID, SpeciesID& maxSpeciesID)
{   // species (input) globally anagenetically speciates
        // -> whole archipelago population diverges from mainland sp
        // can only happen to immigrant species
    vector<int> onWhichIslands = findIsl(speciesID);  // vector with islandIDs
                        // (position in mArchipel) where species is present
    if (onWhichIslands.size() < 2)
        assert(!"Drawn species is present on less than 2 islands. "
                          " Something's wrong.. (global anagenesis)\n");
    // daughter species
    const Species sp = mArchipel[onWhichIslands[0]].findSpecies(speciesID);
    const double birthT = sp.readBirth();
    Species spNew(birthT, speciesID, maxSpeciesID.createNewSpeciesID(), 'A');
    for (auto& iIsl : onWhichIslands) {
        mArchipel[iIsl].goExtinct(speciesID);
        mArchipel[iIsl].addSpecies(spNew);
    }
}

void Archipelago::goGlobalExtinct(const SpeciesID& speciesID)
{   // one species (input) goes exinct on all islands it inhabits
    vector<int> onWhichIslands = findIsl(speciesID);
    if (onWhichIslands.size() < 2)
        assert(!"Drawn species is present on less than 2 islands. "
                "Something's wrong.. (global extinction)\n");
    for (auto& Isl : onWhichIslands) {
        mArchipel[Isl].goExtinct(speciesID);
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
        switch(nextEvent) {
            case event_type::global_cladogenesis:
                assert(static_cast<int>(nextEvent) == 5);
                speciateGlobalClado(speciesID, prng, time, maxSpeciesID);
                break;
            case event_type::global_anagenesis:
                assert(static_cast<int>(nextEvent) == 6);
                speciateGlobalAna(speciesID, maxSpeciesID);
                break;
            case event_type::global_extinction:
                assert(static_cast<int>(nextEvent) == 7);
                goGlobalExtinct(speciesID);
                break;
            default:
                assert(!"Event is not global.\n");  //!OCLINT
                break;
        }
    }
    else if (is_local(nextEvent)) {
        // sample island:
        const int n_islands = getNIslands();
        vector<double> sumRatesPerIsland(n_islands, 0);
        for (int i = 0; i < n_islands; ++i) {
            sumRatesPerIsland[i] += mArchipel[i].extractSumOfRates();
        }
        const int isl = drawDisEvent(sumRatesPerIsland, prng);
        // sample species:
        vector<SpeciesID> aliveSpecies = mArchipel[isl].getSpeciesIDs();
        SpeciesID speciesID = drawUniEvent(aliveSpecies, prng);
        switch(nextEvent)
        {
            case event_type::local_immigration:
            {
                assert(static_cast<int>(nextEvent) == 0);
                mArchipel[isl].immigrate(speciesID, time);
                break;
            }
            case event_type::local_migration:
            {
                assert(static_cast<int>(nextEvent) == 1);
                // save the logarithmic growth terms for all islands per island
                vector<double> vLogs(n_islands, 0);
                for (int j = 0; j < n_islands; ++j) {
                    vLogs[j] = mArchipel[j].returnLogGrowth();
                }
                int destinationIsland = mArchipel[isl].drawMigDestinationIsland(
                        isl, vLogs, initialMigrationRate, prng);
                    // output: position of island in mArchipel where the species
                        // migrates to
                Species newSpecies = mArchipel[isl].findSpecies(speciesID);
                mArchipel[destinationIsland].migrate(newSpecies, time);
                break;
            }
            case event_type::local_cladogenesis:
                assert(static_cast<int>(nextEvent) == 2);
                mArchipel[isl].speciateClado(speciesID, time, maxSpeciesID);
                break;
            case event_type::local_anagenesis:
                assert(static_cast<int>(nextEvent) == 3);
            mArchipel[isl].speciateAna(speciesID, maxSpeciesID);
                break;
            case event_type::local_extinction:
                assert(static_cast<int>(nextEvent) == 4);
            mArchipel[isl].goExtinct(speciesID);
                break;
            default:
                assert(!"Event is not local.\n");  //!OCLINT
                break;
        }
    }
    else
        assert(!"vHappening.size() is neither 2 nor 3. Something is wrong.\n"); //!OCLINT

    // update archi-level extant species vector
//    updateAliveSpec();
}

void Archipelago::addArchi(const Archipelago &newArchi)
{   // adds one archipelago data frame (mArchipel) to this one
    // important for putting together the 1-coloniser-archipelagos
        // -> means, there are no duplicates ?!

    const vector<Island> &addArch = newArchi.returnArchi();
    if(addArch.size() != mArchipel.size())
        throw logic_error("You cannot merge archipelagos of different size.\n");

    // consolidate single islands together
    const int n_islands = getNIslands();
    for (int i = 0; i < n_islands; ++i) {
        if (!addArch[i].returnIsland().empty()) { // ### CAUTION ### : does this make sense ??
            mArchipel[i].consolidateIslands(addArch[i]);
        }
    }
}

vector<Species> Archipelago::makeArchiTo1Island() const
{   // put islands of archipelago together in one Island vector,
        // delete duplicates, sort by time of birth/immigration

    vector<Species> archiToIsland;
    // add all island species vectors together
    for(auto& island : mArchipel) {
        if (!island.returnIsland().empty()) {
            const vector<Species>& vTmp = island.returnIsland();
            archiToIsland.reserve(archiToIsland.size() + vTmp.size());
            archiToIsland.insert(archiToIsland.begin(),
                    vTmp.begin(), vTmp.end());
        }
    }

    // delete duplicates; ### CAUTION ### : what birth time ?!
    for (int j = 0; j < static_cast<int>(archiToIsland.size()); ++j) {
        for (int k = j + 1; k < static_cast<int>(archiToIsland.size()); ++k)
            if (archiToIsland[j].readSpID() ==
                    archiToIsland[k].readSpID()) { // TODO
                // ExtinctTime: take the extant one, or the later extinction
                // BirthTime: take the oldest birth time (initial colonisation) or
                    // the latest re-immigration time.. ### CAUTION ### : How??
                archiToIsland[k] = archiToIsland.back();
                archiToIsland.pop_back();
                --k;
            }
    }
    // sort by birth time TODO
    const int newVecSize = static_cast<int>(archiToIsland.size());
    for (int l = 0; l < newVecSize-1; ++l) {
        if(archiToIsland[l].readBirth() < archiToIsland.back().readBirth()) {
            const Species tmpSp = archiToIsland.back();
            archiToIsland.back() = archiToIsland[l];
            archiToIsland[l] = tmpSp;
        }
    }
    return archiToIsland;
}

void Archipelago::printArchi()
{
    int i = 0;
    for (auto& z : mArchipel) {
        ++i;
        cout << "Island " << i << '\n';
        cout << "BirT" << '\t' << "ParID" << '\t'
             << "SpID" << '\t' << "Stac" << '\n';
        z.printIsland();
    }
}

int Archipelago::getNSpecies()
{
    const int n_aliveSpecies = static_cast<int>(getSpeciesIDs().size());
    return n_aliveSpecies;
}

vector<SpeciesID> Archipelago::getSpeciesIDs()
{
    vector<SpeciesID> aliveSpecies;
    for(auto& island : mArchipel) {
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
    for(auto& island : mArchipel) {
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
    return static_cast<int>(mArchipel.size());
}
