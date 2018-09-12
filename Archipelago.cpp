//
// Created by Bastophiles on 06.09.2018.
//


#include "Archipelago.h"

using namespace std;

Archipelago::Archipelago(const int &nI, const int &AK)  // constructor of archipelago based on number
{                                                            // of islands and archipelago-wide K
    assert(nI >= 0); assert(AK >= 0);
    if(AK % nI != 0)
        throw logic_error("The carrying capacity needs to be dividable"
                          " without remainder by the number of islands.");
    // ### CAUTION ### : let the archipelago-wide K be a double ??
    if (nI == 0 || AK == 0)
        throw string("You are creating an archipelago without islands or a carrying capacity of zero.");
    mvArchipel = vector<Island>((unsigned) nI, Island(AK / nI));
    mAK = AK;
}

/*
void Archipelago::updateAliveSpec()
{   // update archipelago-wide extant species vector
    vector<int> vArchiTmpVec;
    for (auto& i : mvArchipel) {
        if(!i.returnIslSpecAlive().empty()) {
            const vector<int>& vIslTmpVec = i.returnIslSpecAlive();
            vArchiTmpVec.reserve(vArchiTmpVec.size() + vIslTmpVec.size());
            vArchiTmpVec.insert(vArchiTmpVec.end(), vIslTmpVec.begin(), vIslTmpVec.end());
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

vector<int> Archipelago::findIsl(const int &ID) const    // find the island(s) where certain species (input) is within archipelago
                                        // returns vector with island IDs (spot in mvArchipel vector)
// ### CAUTION ### : maybe vector<pairs<int, int> > with island ID and position ??
{
    vector<int> vLocations;
    const int n_islands = mvArchipel.size();
    for (int i = 0; i < n_islands; ++i) {
        if (mvArchipel[i].findPos(ID) >= 0)
            vLocations.push_back(i);
    }
    return vLocations;
}

vector<double> Archipelago::calculateAllRates(const vector<double> &vIniPars, const int &iM, const int &iNumIsl)
{   // the rates for each event on each island are calculated and saved within the island-class
    // order of output: first -> global, second -> local
    // order of parameters (input): gam_i (0), gam_m (1), lamb_cl (2), lamb_al (3), mu_l (4), lamb_cg (5), lamb_ag (6), mu_g (7)

    assert(vIniPars.size() == 8);
    // assign initial parameter values to variables
    const double dIniLamb_Cg = vIniPars[5];
    const double dIniLamb_Ag = vIniPars[6];
    const double dIniMu_g = vIniPars[7];
    const vector<double> vIslPars = { vIniPars[0], vIniPars[1], vIniPars[2], vIniPars[3], vIniPars[4] };
    // order of island parameter vector: gam_i, gam_m, lamb_cl, lamb_al, mu_l
        // for giving it to the island rates calculation function

    // count extant species on archipelago, for global rates
    int iAliveSpecies = getNSpeciesAlive();

    // count immigrant species; only ones that can undergo anagenesis
    int iImmiAlive = 0;
    for (auto& isl : mvArchipel) {
        vector<Species> tmpIsland = isl.returnIsland();
        for (auto& species : tmpIsland)
            if (species.readSpID() <= iM)
                ++iImmiAlive;
    }
    // calculate global rates
    const double dLamb_Cg = max(0.0, dIniLamb_Cg * iAliveSpecies * (1 - (static_cast<double>(iAliveSpecies)/mAK)));   // global cladogenesis
    const double dLamb_Ag = max(0.0, dIniLamb_Ag * iImmiAlive);    // global anagenesis
    const double dMu_g = max(0.0, dIniMu_g * iAliveSpecies);  // global extinction

    const double dSumGlobal = dLamb_Cg + dLamb_Ag + dMu_g;
    mvGlobalRates = { dLamb_Cg, dLamb_Ag, dMu_g };

    // calculate local rates
    // logistic growth term for migration
    double dSumLogGrowth = 0;
    for (auto &j : mvArchipel) {
        dSumLogGrowth += j.returnLogGrowth(); // sums the logistic growth terms of all islands together
    }
    double dSumLocal = 0;
    for (auto &i : mvArchipel) {
        double dThisMigLogGrowth = dSumLogGrowth - i.returnLogGrowth(); // sum of log-growth of all except THIS island
        dSumLocal += i.calculateIslRates(vIslPars, iM, iNumIsl, dThisMigLogGrowth);
    }

    vector<double> pSums = { dSumGlobal, dSumLocal };

    return pSums;
}

vector<int> Archipelago::sampleNextEvent(const vector<double> &pLocalGlobalRates, mt19937_64 prng, const int &iM)
{   // which event, where and to whom will happen next ('when' is calculated in ArchiDAISIE_core)
    // -> all the stochastics happen here

    vector<int> vHappening; // vector for return (global: {event, species} (local: + island)}
    // local or global
    const int iScale = drawDisEvent(pLocalGlobalRates, prng);
    if (iScale == 0) {  // if global event:
        // draw event (plus 5 to align it with all events: 0-4 -> local, 5-7 -> global: 5 -> clado, 6 -> ana, 7 -> extinc)
        const int iEvent = drawDisEvent(mvGlobalRates, prng) + 5;
        // draw species
        vector<int> aliveSpecies;
        for(auto& island : mvArchipel) {
            vector<Species> tmpIsland = island.returnIsland();
            for(auto& species : tmpIsland)
            if (species.isExtant())
                aliveSpecies.push_back(species.readSpID());
        }
        const int pos = drawUniEvent(0, static_cast<int>(aliveSpecies.size() - 1), prng);
        const int iSpecID = aliveSpecies[pos];
        // initialise vector with event and species
        vHappening = { iEvent, iSpecID };
    }
    else {  // -> iScale=1 => if local event
        // draw island
        vector<double> vSumIslRates(mvArchipel.size(),0);
        const int n_islands = mvArchipel.size();
        for (int i = 0; i < n_islands; ++i) {
            vSumIslRates[i] += mvArchipel[i].extractSumIslRate();
        }
        const int iIsl = drawDisEvent(vSumIslRates, prng);
        // initialise vector with event, species and island
        vHappening = mvArchipel[iIsl].sampleLocalEvent(prng, iM);
        vHappening.push_back(iIsl);
    }
    return vHappening;
}

// global updates:

void Archipelago::speciateGlobalClado(const int& speciesID, mt19937_64 prng, double time, SpeciesID& maxSpeciesID)
{   // species (input) globally cladogenetically speciates -> archipelago population splits into two new species
    vector<int> vOnIslands = findIsl(speciesID);  // vector with islandIDs (position in mvArchipel) where species is present
    if (vOnIslands.empty())
        throw logic_error("Drawn species is not present on any island. Something's wrong.. (global cladogenesis)\n");

    // two daughter species
    const int iOnePos = mvArchipel[vOnIslands[0]].findPos(speciesID);
    const double dBirthT = mvArchipel[vOnIslands[0]].returnSpecies(iOnePos).readBirth();
    Species spNew1(dBirthT, speciesID, maxSpeciesID.createNewSpeciesID());
    Species spNew2(time, speciesID, maxSpeciesID.createNewSpeciesID());

    // draw where to split the archipelago: 0 to i-1 -> split after the island number drawn
    int iNumIsl = static_cast<int>(mvArchipel.size());
    const int iSplit = drawUniEvent(0, iNumIsl-1, prng);

    // update data frame
    for (auto& iIsl : vOnIslands) {
        if (iIsl <= iSplit) {
            mvArchipel[iIsl].goExtinct(speciesID, time);
            mvArchipel[iIsl].addSpecies(spNew1);
        }
        else {
            mvArchipel[iIsl].goExtinct(speciesID, time);
            mvArchipel[iIsl].addSpecies(spNew2);
        }
    }
}

void Archipelago::speciateGlobalAna(const int& speciesID, double time, SpeciesID& maxSpeciesID)
{   // species (input) globally anagenetically speciates -> whole archipelago population diverges from mainland sp
    // can only happen to immigrant species
    vector<int> vOnIslands = findIsl(speciesID);  // vector with islandIDs (position in mvArchipel) where species is present
    if (vOnIslands.empty())
        throw logic_error("Drawn species is not present on any island. Something's wrong.. (global anagenesis)\n");
    // daughter species
    const int iOnePos = mvArchipel[vOnIslands[0]].findPos(speciesID);
    const double dBirthT = mvArchipel[vOnIslands[0]].returnSpecies(iOnePos).readBirth();
    Species spNew(dBirthT, speciesID, maxSpeciesID.createNewSpeciesID());
    for (auto& iIsl : vOnIslands) {
        mvArchipel[iIsl].goExtinct(speciesID, time);
        mvArchipel[iIsl].addSpecies(spNew);
    }
}

void Archipelago::goGlobalExtinct(const int& iSpecID, double dTime)
{   // one species (input) goes exinct on all islands it inhabits
    vector<int> vOnIslands = findIsl(iSpecID);
    if (vOnIslands.empty())
        throw logic_error("Drawn species is not present on any island. Something's wrong.. (global extinction)\n");
    for (auto& Isl : vOnIslands) {
        mvArchipel[Isl].goExtinct(iSpecID, dTime);
    }
}

void Archipelago::updateArchi(const vector<int>& happening,
        const double& iniMigrationRate,
        mt19937_64 prng,
        double time,
        SpeciesID& maxSpeciesID)
{   // based on the outcome of sampleNextEvent-function it will update the data frame(s)
    // order of input: event [0], species [1], (island [2])
    // order of parameter indexes (Event): gam_i (0), gam_m (1), lamb_cl (2), lamb_al (3), mu_l (4), lamb_cg (5), lamb_ag (6), mu_g (7)
    assert(happening.size() == 2 || happening.size() == 3);

    const event_type iEvent = static_cast<event_type>(happening[0]);
    const int speciesID = happening[1];

    if (happening.size() == 2) {   // -> global
        //assert(iEvent <= 7 && iEvent >= 5);
        assert(is_global(iEvent));
        switch(iEvent) {
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
                throw logic_error("Event is not global, even though .size() == 2.\n");
        }
    }
    else if (happening.size() == 3) {  // -> local
        assert(is_local(iEvent));
        const int isl = happening[2];
        switch(iEvent)
        {
            case event_type::immigration:
            {
                assert(!"TODO");
                mvArchipel[isl].immigrate(speciesID, time);
                break;
            }
            case event_type::local_migration:
            {
                // ### CAUTION ### : this also includes the island the species is on, right?? -> take care of that !!
                const int n_islands = mvArchipel.size();
                vector<double> vLogs(n_islands);    // minus the island it's on (= isl)
                for (int j = 0; j < n_islands; ++j) {
                    vLogs[j] = mvArchipel[j].returnLogGrowth();
                }
                int destinationIsland = mvArchipel[isl].drawMigDestinationIsland(speciesID,
                        vLogs, iniMigrationRate, prng);
                    // output: position of island in mvArchipel equals island ID
                const Species newSpecies = mvArchipel[isl].findSpecies(speciesID);
                const int pos = mvArchipel[destinationIsland].findPos(speciesID);
                if (pos == -1)       // if first migration: add species to island
                                     // else (if remigration): nothing happens
                    mvArchipel[destinationIsland].addSpecies(newSpecies);   // species (iSpec) migrates
                            // from the original event island to drawn island of destination
                break;
            }
            case event_type::local_cladogenesis:
                mvArchipel[isl].speciateClado(speciesID, time, maxSpeciesID);
                break;
            case event_type::local_anagenesis:
                mvArchipel[isl].speciateAna(speciesID, time, maxSpeciesID);
                break;
            case event_type::local_extinction:
                mvArchipel[isl].goExtinct(speciesID, time);
                break;
            default:
                throw logic_error("Event is not local, even though .size() == 3.\n");
        }
    }
    else
        throw runtime_error("vHappening.size() is neither 2 nor 3. Something is wrong.\n");

    // update archi-level extant species vector
//    updateAliveSpec();
}

void Archipelago::addArchi(const Archipelago &aNewArchi)
{   // adds one archipelago data frame (mvArchipel) to this one (-> member function)
    // important for putting together the 1-coloniser-archipelagos -> means, there are no duplicates ?!
    const vector<Island> &vAddArch = aNewArchi.returnArchi();

    // consolidate single islands together
    const int n_islands = mvArchipel.size();
    for (int i = 0; i < n_islands; ++i) {
        assert(i >= 0);
        assert(i < static_cast<int>(mvArchipel.size()));

        if (!vAddArch[i].returnIsland().empty())    // ### CAUTION ### : does this make sense ??
            assert(i >= 0);
            assert(i < static_cast<int>(mvArchipel.size()));
            mvArchipel[i].addIsland(vAddArch[i]);
    }
}

vector<Species> Archipelago::aggregateArchi() const
{   // put islands of archipelago together in one Island vector, delete duplicates sort by time of event (birth or extinction)
    vector<Species> vAggregatedArchi;

    // add all island species vectors together
    for(auto& i : mvArchipel) {
        const vector<Species>& vTmp = i.returnIsland();
        vAggregatedArchi.reserve(vAggregatedArchi.size() + vTmp.size());
        vAggregatedArchi.insert(vAggregatedArchi.begin(), vTmp.begin(), vTmp.end());
    }

    // delete duplicates; ### CAUTION ### : what birth time ?!
    for (int j = 0; j < static_cast<int>(vAggregatedArchi.size()); ++j) {
        for (int k = j + 1; k < static_cast<int>(vAggregatedArchi.size()); ++k)
            if (vAggregatedArchi[j].readSpID() == vAggregatedArchi[k].readSpID()) {
                // take the oldest birth time (initial colonisation) or the latest re-immigration time.. ### CAUTION ### : How??

                vAggregatedArchi[k] = vAggregatedArchi.back();
                vAggregatedArchi.pop_back();
                --k;
            }
    }
    // sort by birth time
    return vAggregatedArchi;

}

void Archipelago::printArchi()
{
    int i = 0;
    for (auto& z : mvArchipel) {
        ++i;
        cout << "Island " << i << endl;
        cout << "BirthT" << '\t' << "ParentID" << '\t' << "SpecID" << '\t' << "ExtincT" << endl;
        z.printIsland();
    }
}

int Archipelago::getNSpeciesAlive()
{
    int aliveSpeciesCounter = 0;
    for (auto& isl : mvArchipel)
        aliveSpeciesCounter += isl.getNSpeciesAlive();
    return aliveSpeciesCounter;
}
