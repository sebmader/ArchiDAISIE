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
    // ### CAUTION ### : let the archipelago-wide K be a floating point ??
    if (n_islands == 0 || archiCarryingCap == 0)
        throw string("You are creating an archipelago without islands"
                     " and/or with a carrying capacity of zero.\n");
    mArchipel = vector<Island>((unsigned) n_islands,
            Island(archiCarryingCap / n_islands));
    mArchiK = archiCarryingCap;
}

/*
void Archipelago::updateAliveSpec()
{   // update archipelago-wide extant species vector
    vector<int> vArchiTmpVec;
    for (auto& i : mArchipel) {
        if(!i.returnIslSpecAlive().empty()) {
            const vector<int>& vIslTmpVec = i.returnIslSpecAlive();
            vArchiTmpVec.reserve(vArchiTmpVec.size() + vIslTmpVec.size());
            vArchiTmpVec.insert(vArchiTmpVec.end(), vIslTmpVec.begin(),
                vIslTmpVec.end());
        }
    }
    // delete duplicates
    for (size_t i = 0; i < vArchiTmpVec.size(); ++i) {
        for (size_t j = i + 1; j < vArchiTmpVec.size(); ++j)
            if (vArchiTmpVec[i] == vArchiTmpVec[j]) {
                vArchiTmpVec[j] = vArchiTmpVec.back();
                vArchiTmpVec.pop_back();
                --j;
            }
    }
}
*/

vector<int> Archipelago::findIsl(const int &ID) const
// find the island(s) where certain species (input) is within archipelago
// returns vector with island IDs (spots in mArchipel vector)
// ### CAUTION ### : maybe vector<pairs<int, int> > with island ID and position?
{
    vector<int> vLocations;
    const int n_islands = mArchipel.size();
    for (int i = 0; i < n_islands; ++i) {
        if (mArchipel[i].findPos(ID) >= 0)
            vLocations.push_back(i);
    }
    return vLocations;
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

    // count extant species on archipelago, for global rates
    int n_speciesAlive = getNSpeciesAlive();

    // count immigrant species; only ones that can undergo anagenesis
    int n_immigrantsAlive = 0;
    for (auto& isl : mArchipel) {
        vector<Species> tmpIsland = isl.returnIsland();
        for (auto& species : tmpIsland)
            if (species.readSpID() <= n_mainlandSpecies)
                ++n_immigrantsAlive;
    }
    // calculate global rates
    // global cladogenesis:
    const double globalCladoRate = max(0.0, initialGlobalClado * n_speciesAlive *
            (1 - (static_cast<double>(n_speciesAlive)/mArchiK)));
    // global anagenesis:
    const double globalAnaRate = max(0.0, initialGlobalAna * n_immigrantsAlive);
    // global extinction:
    const double globalExtinctRate = max(0.0, initialGlobalExtinct * n_speciesAlive);

    const double sumGlobal = globalCladoRate + globalAnaRate + globalExtinctRate;
    mGlobalRates = { globalCladoRate, globalAnaRate, globalExtinctRate };

    // calculate local rates
    // logistic growth term for migration
    double sumLogGrowth = 0.0;
    for (auto &j : mArchipel) {
        sumLogGrowth += j.returnLogGrowth(); // sums the
            // logistic growth terms of all islands together
    }
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

vector<int> Archipelago::sampleNextEvent(const vector<double> &sumsGloLoc,
        mt19937_64 prng,
        const int &n_mainlandSpecies)
{   // which event, where and to whom will happen next
        // ('when' is calculated in ArchiDAISIE_core)
    // -> all the stochastics happen here

    vector<int> happening; // vector for return
        // global: {event, species}; local: {event, species, island}

    // local or global:
    const int scale = drawDisEvent(sumsGloLoc, prng);
    if (scale == 0) {  // if global event:

        // draw event (plus 5 to align it with all events:
            // 0-4 -> local, 5-7 -> global: 5 -> clado, 6 -> ana, 7 -> extinc)
        const int event = drawDisEvent(mGlobalRates, prng) + 5;

        // draw species
        vector<int> aliveSpecies = getIDsSpeciesAlive();
        const int maxPosAlive = static_cast<int>(aliveSpecies.size() - 1);
        const int pos = drawUniEvent(0,
                maxPosAlive, prng);
        const int speciesID = aliveSpecies[pos];

        // initialise vector with event and species
        happening = { event, speciesID };
    }
    else {  // -> scale=1 => if local event
        // draw island
        const int n_islands = mArchipel.size();

        vector<double> sumRatesPerIsland(n_islands,0);
        for (int i = 0; i < n_islands; ++i) {
            sumRatesPerIsland[i] += mArchipel[i].extractSumOfRates();
        }
        const int island = drawDisEvent(sumRatesPerIsland, prng);

        // initialise vector with event, species and island
        happening = mArchipel[island].sampleLocalEvent(prng,
                n_mainlandSpecies);
        happening.push_back(island);
    }
    return happening;
}

// global updates:

void Archipelago::speciateGlobalClado(const int& speciesID,
        mt19937_64 prng,
        double time,
        SpeciesID& maxSpeciesID)
    // species (input) globally cladogenetically speciates
{   // -> archipelago population splits into two new species
    vector<int> onWhichIslands = findIsl(speciesID);  // vector with islandIDs
        // (position in mArchipel) where species is present
    if (onWhichIslands.empty())
        assert(!"Drawn species is not present on any island. "
                          "Something's wrong.. (global cladogenesis)\n");

    // two daughter species
    const Species sp = mArchipel[onWhichIslands[0]].findSpecies(speciesID);
    const double birthT = sp.readBirth();
    Species spNew1(birthT, speciesID, maxSpeciesID.createNewSpeciesID());
    Species spNew2(time, speciesID, maxSpeciesID.createNewSpeciesID());

    // draw where to split the archipelago:
        // 0 to i-2 -> split after the island number drawn
    int n_islands = static_cast<int>(mArchipel.size());
    const int split = drawUniEvent(0, n_islands-2, prng);

    // update data frame
    for (auto& isl : onWhichIslands) {
        mArchipel[isl].goExtinct(speciesID, time);
        if (isl <= split)
            mArchipel[isl].addSpecies(spNew1);
        else
            mArchipel[isl].addSpecies(spNew2);
    }
}

void Archipelago::speciateGlobalAna(const int& speciesID,
        double time,
        SpeciesID& maxSpeciesID)
{   // species (input) globally anagenetically speciates
        // -> whole archipelago population diverges from mainland sp
        // can only happen to immigrant species
    vector<int> onWhichIslands = findIsl(speciesID);  // vector with islandIDs
                        // (position in mArchipel) where species is present
    if (onWhichIslands.empty())
        assert(!"Drawn species is not present on any island."
                          " Something's wrong.. (global anagenesis)\n");
    // daughter species
    const Species sp = mArchipel[onWhichIslands[0]].findSpecies(speciesID);
    const double birthT = sp.readBirth();
    Species spNew(birthT, speciesID, maxSpeciesID.createNewSpeciesID());
    for (auto& iIsl : onWhichIslands) {
        mArchipel[iIsl].goExtinct(speciesID, time);
        mArchipel[iIsl].addSpecies(spNew);
    }
}

void Archipelago::goGlobalExtinct(const int& speciesID, double time)
{   // one species (input) goes exinct on all islands it inhabits
    vector<int> onWhichIslands = findIsl(speciesID);
    if (onWhichIslands.empty())
        assert(!"Drawn species is not present on any island. "
                "Something's wrong.. (global extinction)\n");
    for (auto& Isl : onWhichIslands) {
        mArchipel[Isl].goExtinct(speciesID, time);
    }
}

void Archipelago::doNextEvent(const vector<int>& happening,
        const double& initialMigrationRate,
        mt19937_64 prng,
        double time,
        SpeciesID& maxSpeciesID)
{   // based on the outcome of sampleNextEvent function
        // it will update the data frame(s)
    // order of input: event [0], species [1], (island [2])
    // order of parameter indexes (Event):
        // immigration (0), migration (1), localClado(2), localAna (3),
        // localExtinct (4), globalClado(5), GlobalAna(6), GlobalExtinct(7)
    assert(happening.size() == 2 || happening.size() == 3);

    const auto event = static_cast<event_type>(happening[0]);
    const int speciesID = happening[1];

    if (happening.size() == 2) {   // -> global
        assert(is_global(event));
        switch(event) {
            case event_type::global_cladogenesis:
                speciateGlobalClado(speciesID, prng, time, maxSpeciesID);
                break;
            case event_type::global_anagenesis:
                speciateGlobalAna(speciesID, time, maxSpeciesID);
                break;
            case event_type::global_extinction:
                goGlobalExtinct(speciesID, time);
                break;
            default:
                assert(!"Event is not global, even though .size() == 2.\n");
                break;
        }
    }
    else if (happening.size() == 3) {  // -> local
        assert(is_local(event));
        const int isl = happening[2];
        switch(event)
        {
            case event_type::immigration:
            {
                assert(!"TODO");
                mArchipel[isl].immigrate(speciesID, time);
                break;
            }
            case event_type::local_migration:
            {
                const int n_islands = mArchipel.size();
                // save the logarithmic growth terms for all islands per island
                vector<double> vLogs(n_islands);
                for (int j = 0; j < n_islands; ++j) {
                    vLogs[j] = mArchipel[j].returnLogGrowth();
                }
                int destinationIsland = mArchipel[isl].drawMigDestinationIsland(
                        isl, vLogs, initialMigrationRate, prng);
                    // output: position of island in mArchipel whereto the species
                        // migrates
                const Species newSpecies = mArchipel[isl].findSpecies(speciesID);
                const int pos = mArchipel[destinationIsland].findPos(speciesID);
                if (pos == -1)  // if first migration: add species to island
                                // else (if remigration): nothing happens
                    mArchipel[destinationIsland].addSpecies(newSpecies);   // species migrates
                            // from the original event island to drawn island of destination
                break;
            }
            case event_type::local_cladogenesis:
                mArchipel[isl].speciateClado(speciesID, time, maxSpeciesID);
                break;
            case event_type::local_anagenesis:
                mArchipel[isl].speciateAna(speciesID, time, maxSpeciesID);
                break;
            case event_type::local_extinction:
                mArchipel[isl].goExtinct(speciesID, time);
                break;
            default:
                assert(!"Event is not local, even though .size() == 3.\n");
                break;
        }
    }
    else
        assert(!"vHappening.size() is neither 2 nor 3. Something is wrong.\n");

    // update archi-level extant species vector
//    updateAliveSpec();
}

void Archipelago::addArchi(const Archipelago &newArchi)
{   // adds one archipelago data frame (mArchipel) to this one
    // important for putting together the 1-coloniser-archipelagos
        // -> means, there are no duplicates ?!

    const vector<Island> &addArch = newArchi.returnArchi();

    // consolidate single islands together
    const int n_islands = static_cast<int>(mArchipel.size());
    for (int i = 0; i < n_islands; ++i) {
        assert(i >= 0);
        assert(i < n_islands);

        if (!addArch[i].returnIsland().empty()) { // ### CAUTION ### : does this make sense ??
            mArchipel[i].addIsland(addArch[i]);
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
    const int vecSize = static_cast<int>(archiToIsland.size());
    for (int j = 0; j < vecSize; ++j) {
        for (int k = j + 1; k < vecSize; ++k)
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
        cout << "BirthT" << '\t' << "ParentID" << '\t'
             << "SpecID" << '\t' << "ExtincT" << '\n';
        z.printIsland();
        cout << '\n';
    }
}

int Archipelago::getNSpeciesAlive()
{
    int aliveSpeciesCounter = 0;
    for (auto& isl : mArchipel)
        aliveSpeciesCounter += isl.getNSpeciesAlive();
    return aliveSpeciesCounter;
}

vector<int> Archipelago::getIDsSpeciesAlive()
{
    vector<int> aliveSpecies;
    for(auto& island : mArchipel) {
        vector<int> tmpIDs = island.getIDsSpeciesAlive();
        aliveSpecies.reserve(aliveSpecies.size() + tmpIDs.size());
        aliveSpecies.insert(aliveSpecies.end(), tmpIDs.begin(), tmpIDs.end());
    }
    return aliveSpecies;
}
