/*  ArchiDAISIE - Simulation of species assembly on an
 *                archipelago from emergence to present
 *  written by Sebastian Mader (S3503704) - 10-07-2018
*/

// ------------ INCLUDED LIBRARIES ------------ //

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <string>
#include <vector>
#include <fstream>
#include <utility>
#include <random>
#include <chrono>
#include <exception>

using namespace std;

// ------------ GLOBAL PARAMETERS ------------ //

    // do I even need globals?
    // should be provided as function parameters, right?

    double dTime;

// ------------ FUNCTION DEFINITIONS ------------ //

// --- distribution functions ---
inline int drawDisEvent(const vector<double> &vecRates, mt19937_64 &prng) {
//draw a discrete distribution event
    discrete_distribution<int> drawEvent(vecRates.begin(), vecRates.end());
    return drawEvent(prng);
}

inline int drawUniEvent(const int &botBoundary, const int &topBoundary, mt19937_64 &prng) {
//draw a uniform distribution event
    // using "topBoundary" (= size of vector - 1 = position of last element of vector)
    // because this function works with different types of vectors: vector<int> & vector<pair<int,int> >
    uniform_int_distribution<int> drawEvent(0, topBoundary);
    return drawEvent(prng);
}


// ------------ CLASS DEFINITIONS ------------ //

#include "Species.h"


class Archipelago {         // class for the whole archipelago
public:
    class Island {        // class for ONE island within archipelago
    public:
        explicit Island(const int k) : mIK{k} {assert(k >= 0);} // island constructor based on island-wide K

        // int sizeIsl() const {return static_cast<int>(mvIsland.size());}  // returns size of island vector
        int specAlive() const {return static_cast<int>(mvIslSpecAlive.size());} // returns number of species alive
        double returnLogGrowth() { return 1 - static_cast<double>(mvIslSpecAlive.size()) / mIK;}   // returns the logistic growth term (1-n/K)

        int findPos(const int &ID) const;    // find the position of certain species (input) in IslandPhylo vector

        const int createNewID();       // returns new species ID and maxID += 1

        double calculateIslRates(const vector<double> &, const int &, const int &, const double &);
                                            // initialise/calculate rates and store them in EventRates vector
                                            // gam_i, gam_m, lamb_cl, lamb_al, mu_l
                                            // per island -> doesn't include global rates !!!
                                            const double extractSumIslRate() const noexcept;      // return the per-island rates vector

        vector<int> sampleLocalEvent(mt19937_64, const int&);   // in case a local event is drawn, sample island, event and species
                                            // it happens to

        void immigrate(const int&, const double& BirthT);                   // mainland species immigrates to that island
        int migrate(const int &, vector<double> &, const double &, mt19937_64);                     // island species migrates to other island
        void speciateClado(const int &);               // island species cladogenetically speciates
        void speciateAna(const int &);                 // island species anagenetically speciates
        void goExtinct(const int &);                   // island species goes extinct

        Species returnSpecies(const int &iPos) { return mvIsland[iPos]; }   // returns specific species from species vector
        void pushbackSp(const Species &spNew) { mvIsland.push_back(spNew); }    // adds new species to species vector

        void addIsland(const Island &);     // add island to THIS islandd

        void printIsland();                 // prints island vector of species to the screen

        const vector<Species>& returnIsland() const { return mvIsland; }    // return island vector of species
        const vector<int>& returnIslSpecAlive() const { return mvIslSpecAlive; }  // return extant species vector
    private:
        vector<Species> mvIsland;           // phylogeny vector of species on island
        vector<int> mvIslSpecAlive;            // vector of species' IDs of alive species (.size() = n species on island)
        vector<double> mvLocalRates;        // vector of rates for events PER ISLAND (5 rates: gamI, gamM, lambC, lambA, mu)
        int mIK; // Carrying capacity (should be const one day)
                 // for now: mIk = mAK / iNumIslands
    };

    Archipelago(const int &, const int &);  // constructor of archipelago based on number of islands and
                                            // archipelago-wide K
    void updateAliveSpec();             // updating the archipelago-wide vector of extant species

    vector<int> findIsl(const int &) const;    // find the island(s) where certain species (input) is within archipelago

    vector<double> calculateAllRates(const vector<double> &, const int &iM, const int &iNumIsl);    // calculate per-island rates and global rates
                                            // and save them in LocalRates and GlobalRates vector, resp.
                                            // Also, output of sum of both global (.first) and local (.second) rates
    vector<int> sampleNextEvent(const vector<double> &, mt19937_64, const int &);   // draw next event; output -> {event(0-7), species(ID)(,island(0-i))}
                                        // if global event -> vector.size() = 2, if local -> size = 3
    const int createNewID();       // returns new species ID and maxID += 1

    void speciateGlobalClado(const int &, mt19937_64);         // island species cladogenetically speciates over all islands
                                        // (one population gets replaced by two new species populations)
                                        // with random separation of archipelago into two groups/populations
    void speciateGlobalAna(const int &);           // island species collectively (on all islands) diverges from mainland ancestor
    void goGlobalExtinct(const int &);             // island species goes extinct on all islands it occures on

    void updateArchi(const vector<int> &, const double &, mt19937_64); // switch-statement that calls event functions
                                        // updates the ArchiPhylo vector
                                        // LOCAL events indicated by 3 elements: event ([0]), species ([1]), island ([2])
                                        // GLOBAL events indicated by 2 elements: event ([0]), species ([1])
    void addArchi(const Archipelago&);  // add an ArchiPhylo to this one / consolidate them

    const vector<Species>& aggregateArchi();   // aggregate all islands in ArchiPylo vector as it would be one and return it
    const vector<Island> & returnArchi() const {return mvArchipel;}      // return OneIslaArchiPylo vector (and output to file ?)
                        // ### CAUTION ###: only const& if assigning this output to a variable copies it
    void printArchi();

    static void setMaxID(const int &maxID) {Archipelago::mMaxSpID = maxID;}
    static int& returnMaxID() {return Archipelago::mMaxSpID;}
    static void incrementMaxID() {++Archipelago::mMaxSpID;}

private:
    vector<double> mvGlobalRates;           // vector of rates for global events (3 rates: lambC, lambA, mu)
    vector<Island> mvArchipel;            // vector of island phylogenies of whole archipelago
    vector<int> mvArchSpecAlive;        // vector of IDs of extant species
                        // ### CAUTION ### : how to update? -> at the end of each event / within the event update functions
    // vector<Species> mvOneIslArchipel;    // aggregated vector of phylogeny of whole archipelago (like it'd be one island)
    static int mMaxSpID;    // the currently highest species ID
            // ### CAUTION ### : static! --> one shared member var for all archipelago objects initialised
    int mAK;    // carrying capacity of archipelago
};

int Archipelago::mMaxSpID = 0;  // defining the static int MaxSp (on archipelago!) as 0 for the start

// ------------ Abbreviations ------------ //

using Isl = Archipelago::Island;
using Archi = Archipelago;


// ------------ Class Member Functions ------------ //

// --- Island: ---

void Archipelago::Island::printIsland()
{
    for (auto &i : mvIsland)
        i.printSpec();
}

int Archipelago::Island::findPos(const int &ID) const    // find the position of certain species (input) in species vector
{                                                        // if not present: output = -1
    for (int i = 0; i < mvIsland.size(); ++i)
        if (mvIsland[i].readSpID() == ID)
            return i;
    return -1;
}

const int Archipelago::Island::createNewID()
{
    incrementMaxID();
    return returnMaxID();
}

double Archipelago::Island::calculateIslRates(const vector<double> &vIslPars, const int &iM, const int &iNumIsl,
                                              const double &dThisMLG)
{   // calculates the per-island rates of events, outputs them (for immidiate use)
    // and saves them within the island-class; input -> initial parameters
    // order of island parameter vector: gam_i, gam_m, lamb_cl, lamb_al, mu_l

    // calculate the per-island rates
    const int iN = static_cast<int>(mvIslSpecAlive.size());
    int iNImmi = 0;    // number of immigrant species -> only ones that can undergo anagenesis
    for (int g : mvIslSpecAlive)
        if (g <= 1000)
            ++iNImmi;
    const double dGam_I = max(0.0, (vIslPars[0] * iM * (1 - static_cast<double>(iN) / mIK)) / iNumIsl);  // immigration to specific island
    const double dGam_M = max(0.0, (vIslPars[1] * iN * dThisMLG) / iNumIsl*iNumIsl - iNumIsl); // migration from this island to all other islands
                                    // all other islands: '/iNumIsl', one-way: '/iNumIsl^2 - iNumIsl' -> ### CAUTION ### : really ???
    const double dLamb_Cl = max(0.0, vIslPars[2] * iN * (1 - (static_cast<double>(iN) / mIK))); // local cladogenesis
    const double dLamb_Al = max(0.0, vIslPars[3] * iNImmi);   // local anagenesis
    const double dMu_l = max(0.0, vIslPars[4] * iN);  // local extinction

    // save in rate vector
    mvLocalRates = { dGam_I, dGam_M, dLamb_Cl, dLamb_Al, dMu_l };

    double dSumIslRates = dGam_I + dGam_M + dLamb_Cl + dLamb_Al + dMu_l;

    return dSumIslRates;
}

const double Archipelago::Island::extractSumIslRate() const noexcept
{
    double dSumRates = 0;
    for (double mvLocalRate : mvLocalRates)
        dSumRates += mvLocalRate;
    return dSumRates;
}

vector<int> Archipelago::Island::sampleLocalEvent(mt19937_64 prng, const int &iM)
{   // samples local event and species it happens to on THIS island
    // draw event
    const int iEvent = drawDisEvent(mvLocalRates, prng);
    // draw species
    int iSpecID;
    if (iEvent) {   // if not immigration (1-4) -> SpecID from extant island species
        iSpecID = mvIslSpecAlive[drawUniEvent(0, static_cast<int>(mvIslSpecAlive.size()), prng)];
    }
    else {  // if immigration (0) -> SpecID from mainland species pool
        iSpecID = drawUniEvent(1, iM, prng);
    }

    // return event and specID
    vector<int> vHappening = { iEvent, iSpecID };

    return vHappening;
}

// local updates:
void Archipelago::Island::immigrate(const int& iSpecID, const double& BirthT = dTime)
{   // immigration from the mainland to THIS island
    // check if species is already present on island
    const int iPos = findPos(iSpecID);
    // if mainland sp -> BirthT = Time; else (if island sp) -> BirthT = old BirthT (because this function is used for migration as well !!)
    Species spNew(BirthT, iSpecID, iSpecID);    // immigrant
    if (BirthT != dTime)                        // if not an immigrant:
        spNew = Species(BirthT, iSpecID, createNewID());
    if (iPos == -1) {    // not present
        mvIsland.push_back(spNew);
        mvIslSpecAlive.push_back(iSpecID);
    }
    else {  // if present
        if (!mvIsland[iPos].isExtant()) {   // if extinct
            mvIsland.push_back(spNew);
            mvIslSpecAlive.push_back(iSpecID);
        }
        else  // if extant -> re-immigration ("re-setting the clock" (= BirthT))
            mvIsland[iPos] = spNew;
    }
}

int Archipelago::Island::migrate(const int &iSpecID, vector<double> &vLogs, const double &dIniMigRate, mt19937_64 prng)
{   // migration from THIS island to another; output: island of destination
    // draw island to which species migrates -> initial migration rate as parameter !!
    const unsigned long iNumIsl = vLogs.size();
    vector<double> vMigRates(iNumIsl);
    for (int k = 0; k < iNumIsl; ++k) {
        vMigRates[k] = max(0.0, (dIniMigRate * mvIslSpecAlive.size() * vLogs[k]) / iNumIsl*iNumIsl - iNumIsl);
    }
    const int iDestinationIsl = drawDisEvent(vMigRates, prng);

    return iDestinationIsl;
}

void Archipelago::Island::speciateClado(const int &iSpecID)
{   // island species cladogenetically diverges
    // find species
    const int iPos = findPos(iSpecID);
    if (iPos == -1)
        throw logic_error("Drawn species is not present on island.. Something's wrong (Cladogenesis).\n");
    // 2 new species:
    Species spNew1 = Species(mvIsland[iPos].readBirth(), iSpecID, createNewID());
    Species spNew2 = Species(dTime, iSpecID, createNewID());
    // parent goes extinct
    mvIsland[iPos].goExtinct(dTime);
    mvIsland.push_back(spNew1);
    mvIsland.push_back(spNew2);
}

void Archipelago::Island::speciateAna(const int &iSpecID)
{   // anagenetic speciation: only immigrants can undergo anagenesis !!!
    // find species
    const int iPos = findPos(iSpecID);
    if (iPos == -1)
        throw logic_error("Drawn species is not present on island.. Something's wrong (Anagenesis).\n");
    // new species
    Species spNew = Species(mvIsland[iPos].readBirth(), iSpecID, createNewID());
    // parent goes extinct & daugther gets added to island
    mvIsland[iPos].goExtinct(dTime);
    mvIsland.push_back(spNew);
}

void Archipelago::Island::goExtinct(const int &iSpecID)
{   // species goes extinct
    // find species
    const int iPos = findPos(iSpecID);
    if (iPos == -1)
        throw logic_error("Drawn species is not present on island.. Something's wrong (Extinction).\n");
    // extinction
    mvIsland[iPos].goExtinct(dTime);
}

void Archipelago::Island::addIsland(const Island &islNew)
{   // adds another island to THIS (for aggregating archipelagos)
    // extract data frames from island that's to be added
    const vector<Species>& vSpecNew = islNew.returnIsland();
    const vector<int>& vSpecAliveNew = islNew.returnIslSpecAlive();

    // add species vector to THIS island
    if (!vSpecNew.empty()) {
        // intermediate vector of species
        vector<Species> vTempSpec = mvIsland;
        vTempSpec.reserve(mvIsland.size() + vSpecNew.size());   // ### CAUTION ### : problem is the const type of ALL outputs and functions
                                                                                //  associated with .reserve()
        vTempSpec.insert(vTempSpec.end(), vSpecNew.begin(), vSpecNew.end());
        // mvIsland.reserve(mvIsland.size() + vSpecNew.size());   // preallocate memory
        // mvIsland.insert(mvIsland.end(), vSpecNew.begin(), vSpecNew.end());
    }
    // add alive species to THIS island
    if (!vSpecAliveNew.empty()) {
        vector<int> vTempAlive = mvIslSpecAlive;
        vTempAlive.reserve(mvIslSpecAlive.size() + vSpecAliveNew.size());
        vTempAlive.insert(vTempAlive.end(), vSpecAliveNew.begin(), vSpecAliveNew.end());
        // mvIslSpecAlive.reserve(mvIslSpecAlive.size() + vSpecAliveNew.size());
        // mvIslSpecAlive.insert(mvIslSpecAlive.end(), vSpecAliveNew.begin(), vSpecAliveNew.end());
    }
}

// --- Archipelago: ---

Archipelago::Archipelago(const int &nI, const int &AK)  // constructor of archipelago based on number
{                                                            // of islands and archipelago-wide K
    assert(nI >= 0); assert(AK >= 0); assert(AK % nI == 0);   // ### CAUTION ### : let the archipelago-wide K be a double ??
    if (nI == 0 || AK == 0)
        throw string("You are creating an archipelago without islands or a carrying capacity of zero.");
    mvArchipel = vector<Island>((unsigned) nI, Island(AK / nI));
    mAK = AK;
}

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

vector<int> Archipelago::findIsl(const int &ID) const    // find the island(s) where certain species (input) is within archipelago
                                        // returns vector with island IDs (spot in mvArchipel vector)
// ### CAUTION ### : maybe vector<pairs<int, int> > with island ID and position ??
{
    vector<int> vLocations;
    int i;
    for (i = 0; i < mvArchipel.size(); ++i) {
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
    int iAliveSpecies = static_cast<int>(mvArchSpecAlive.size());

    // count immigrant species; only ones that can undergo anagenesis
    int iImmiAlive = 0;
    for (auto& h : mvArchSpecAlive)
        if (h <= iM)
            ++iImmiAlive;

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
        const int iSpecID = mvArchSpecAlive[drawUniEvent(0, static_cast<int>(mvArchSpecAlive.size() - 1),
                                                         prng)];
        // initialise vector with event and species
        vHappening = { iEvent, iSpecID };
    }
    else {  // -> iScale=1 => if local event
        // draw island
        vector<double> vSumIslRates(mvArchipel.size(),0);
        for (int i = 0; i < mvArchipel.size(); ++i) {
            vSumIslRates[i] += mvArchipel[i].extractSumIslRate();
        }
        const int iIsl = drawDisEvent(vSumIslRates, prng);
        // initialise vector with event, species and island
        vHappening = mvArchipel[iIsl].sampleLocalEvent(prng, iM);
        vHappening.push_back(iIsl);
    }

    return vHappening;
}

const int Archipelago::createNewID()
{
    incrementMaxID();
    return returnMaxID();
}

// global updates:
void Archipelago::goGlobalExtinct(const int &iSpecID)
{   // one species (input) goes exinct on all islands it inhabits
    vector<int> vOnIslands = findIsl(iSpecID);
    if (vOnIslands.empty())
        throw logic_error("Drawn species is not present on any island. Something's wrong.. (global extinction)\n");
    for (auto& Isl : vOnIslands) {
        mvArchipel[Isl].goExtinct(iSpecID);
    }
}

void Archipelago::speciateGlobalAna(const int &iSpecID)
{   // species (input) globally anagenetically speciates -> whole archipelago population diverges from mainland sp
    // can only happen to immigrant species
    vector<int> vOnIslands = findIsl(iSpecID);  // vector with islandIDs (position in mvArchipel) where species is present
    if (vOnIslands.empty())
        throw logic_error("Drawn species is not present on any island. Something's wrong.. (global anagenesis)\n");
    // daughter species
    const int iOnePos = mvArchipel[vOnIslands[0]].findPos(iSpecID);
    const double dBirthT = mvArchipel[vOnIslands[0]].returnSpecies(iOnePos).readBirth();
    Species spNew(dBirthT, iSpecID, createNewID());
    for (auto& iIsl : vOnIslands) {
        mvArchipel[iIsl].goExtinct(iSpecID);
        mvArchipel[iIsl].pushbackSp(spNew);
    }
}

void Archipelago::speciateGlobalClado(const int &iSpecID, mt19937_64 prng)
{   // species (input) globally cladogenetically speciates -> archipelago population splits into two new species
    vector<int> vOnIslands = findIsl(iSpecID);  // vector with islandIDs (position in mvArchipel) where species is present
    if (vOnIslands.empty())
        throw logic_error("Drawn species is not present on any island. Something's wrong.. (global cladogenesis)\n");

    // two daughter species
    const int iOnePos = mvArchipel[vOnIslands[0]].findPos(iSpecID);
    const double dBirthT = mvArchipel[vOnIslands[0]].returnSpecies(iOnePos).readBirth();
    Species spNew1(dBirthT, iSpecID, createNewID());
    Species spNew2(dTime, iSpecID, createNewID());

    // draw where to split the archipelago: 0 to i-1 -> split after the island number drawn
    int iNumIsl = static_cast<int>(mvArchipel.size());
    const int iSplit = drawUniEvent(0, iNumIsl-1, prng);

    // update data frame
    for (auto& iIsl : vOnIslands) {
        if (iIsl <= iSplit) {
            mvArchipel[iIsl].goExtinct(iSpecID);
            mvArchipel[iIsl].pushbackSp(spNew1);
        }
        else {
            mvArchipel[iIsl].goExtinct(iSpecID);
            mvArchipel[iIsl].pushbackSp(spNew2);
        }
    }
}


void Archipelago::updateArchi(const vector<int> &vHappening, const double &dIniMigRate, mt19937_64 prng)
{   // based on the outcome of sampleNextEvent-function it will update the data frame(s)
    // order of input: event [0], species [1], (island [2])
    // order of parameter indexes (Event): gam_i (0), gam_m (1), lamb_cl (2), lamb_al (3), mu_l (4), lamb_cg (5), lamb_ag (6), mu_g (7)

    const int iEvent = vHappening[0];
    const int iSpecID = vHappening[1];

    if (vHappening.size() == 2) {   // -> global
        assert(iEvent <= 7 && iEvent >= 5);
        switch(iEvent) {
            case 5:
                speciateGlobalClado(iSpecID, prng);
                break;
            case 6:
                speciateGlobalAna(iSpecID);
                break;
            case 7:
                goGlobalExtinct(iSpecID);
                break;
            default:
                throw logic_error("Event is not global, even though .size() == 2.\n");
        }
    }
    else if (vHappening.size() == 3) {  // -> local
        assert(iEvent >= 0 && iEvent <= 4);
        const int iIsl = vHappening[2];
        switch(iEvent) {
            case 0:
                mvArchipel[iIsl].immigrate(iSpecID);
                break;
            case 1: {
                vector<double> vLogs(mvArchipel.size());
                for (int j = 0; j < mvArchipel.size(); ++j)
                    vLogs[j] = mvArchipel[j].returnLogGrowth();
                int iDestinationIsl = mvArchipel[iIsl].migrate(iSpecID, vLogs, dIniMigRate, prng);    // output: position of island in mvArchipel
                                                                                                    // equals island ID
                const int iPos = mvArchipel[iIsl].findPos(iSpecID);
                const double dBirthT = mvArchipel[iIsl].returnSpecies(iPos).readBirth();
                mvArchipel[iDestinationIsl].immigrate(iSpecID, dBirthT);   // species (iSpec) "immigrates" from the original event island
                                                                            // to drawn island of destination
                break;
            }
            case 2:
                mvArchipel[iIsl].speciateClado(iSpecID);
                break;
            case 3:
                mvArchipel[iIsl].speciateAna(iSpecID);
                break;
            case 4:
                mvArchipel[iIsl].goExtinct(iSpecID);
                break;
            default:
                throw logic_error("Event is not local, even though .size() == 3.\n");
        }
    }
    else
        throw runtime_error("vHappening.size() is neither 2 nor 3. Something is wrong.\n");

    // update archi-level extant species vector
    updateAliveSpec();
}

void Archipelago::addArchi(const Archipelago &aNewArchi)
{   // adds one archipelago data frame (mvArchipel) to this one (-> member function)
    // important for putting together the 1-coloniser-archipelagos -> means, there are no duplicates ?!
    const vector<Island> &vAddArch = aNewArchi.returnArchi();

    // consolidate single islands together
    for (int i = 0; i < mvArchipel.size(); ++i) {
        if (!vAddArch[i].returnIsland().empty())    // ### CAUTION ### : does this make sense ??
            mvArchipel[i].addIsland(vAddArch[i]);
    }
}

const vector<Species>& Archipelago::aggregateArchi()
{   // put islands of archipelago together in one Island vector, delete duplicates sort by time of event (birth or extinction)
    vector<Species> vAggregatedArchi;

    // add all island species vectors together
    for(auto& i : mvArchipel) {
        const vector<Species>& vTmp = i.returnIsland();
        vAggregatedArchi.reserve(vAggregatedArchi.size() + vTmp.size());
        vAggregatedArchi.insert(vAggregatedArchi.begin(), vTmp.begin(), vTmp.end());
    }

    // delete duplicates; ### CAUTION ### : what birth time ?!
    for (int j = 0; j < vAggregatedArchi.size(); ++j) {
        for (int k = j + 1; k < vAggregatedArchi.size(); ++k)
            if (vAggregatedArchi[j].readSpID() == vAggregatedArchi[k].readSpID()) {
                // take the oldest birth time (initial colonisation) or the latest re-immigration time.. ### CAUTION ### : How??

                vAggregatedArchi[k] = vAggregatedArchi.back();
                vAggregatedArchi.pop_back();
                --k;
            }
    }
    // sort by birth time


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


// ------------ ArchiDAISIE FUNCTIONS ------------ //

Archipelago ArchiDAISIE_core(const double &dAge, const unsigned long int &M, const vector<double> &vIniPars,
                                const int &iAK, const unsigned int &iNumIslands, mt19937_64 &prng)
{
    try {
        // initialise Archipelago data frame and set time and max species ID to initial values
        Archipelago aArchi(iNumIslands, iAK);
        dTime = dAge;

        // initialise a mainland species vector to sample from for immigration
            // PROBABLY: not even needed; you can sample from uniform distribution of 1 to M
            // and number drawn = species ID ..?
/*        vector<int> vMainlandSp(M);
        for (unsigned i = 0; i < vMainlandSp.size(); ++i) {
            vMainlandSp[i] = i + 1;
        }*/


        // start looping through time
        for (;;) {

            // calculate the rates of events
            vector<double> pLAndGRates = aArchi.calculateAllRates(vIniPars, (int) M, iNumIslands);

            // draw time interval to next event
            const double dSumRates = pLAndGRates[0] + pLAndGRates[1];
            if (dSumRates <= 0)
                throw runtime_error("Event rate is zero or below. No event can be drawn.\n");

            // draw time to next event and take it off the actual time
            exponential_distribution<double> waitingTime(dSumRates);
            const double dDT = waitingTime(prng);
            dTime -= dDT;
            if (dTime <= 0)     // when dTime passes the present, stop simulation
                break;

            // sample which event happens
            vector<int> vHappening = aArchi.sampleNextEvent(pLAndGRates, prng, (int) M);

            // update the phylogeny
            aArchi.updateArchi(vHappening, vIniPars[1], prng);
        }
    }
    catch (string &str) {
        std::cerr << "Warning: " << str;
    }
}

vector<vector<Species> > ArchiDAISIE(const double &dAge, const unsigned int iMainSp_n, vector<double> &vIniPars,
                                        const unsigned int &iNumIslands, const unsigned int iReplicates)
{   // ### CAUTION ### : output unclear !!
    try {
        // check given parameters
        if (dAge <= 0)
            throw logic_error("Age has to be higher than zero.");
        if (iMainSp_n <= 0)
            throw logic_error("Simulation needs at least one mainland species.");
        if(iNumIslands <= 0)
            throw logic_error("Simulation needs at least one island.");
        if(vIniPars.size() != 9)
            throw logic_error("Provide 9 parameter values.");
        if(vIniPars[0] <= 0)
            throw logic_error("Rate of colonisation is zero or below. The island cannot be colonised.");

        // declare and seed PRNG with system clock
        mt19937_64 prng;

        chrono::high_resolution_clock::time_point tp =
                chrono::high_resolution_clock::now();
        const unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
        // ### CAUTION: ###   should print/output seed here
        prng.seed(seed);

        // order of parameters (input): gam_i, gam_m, lamb_cl, lamb_al, mu_l, lamb_cg, lamb_ag, mu_g, ArchiK
        const int iAK = static_cast<int>(vIniPars[8]);
        vIniPars.pop_back();

        // initialise main data frame
        vector< vector<Species> > vFinalIslandReplicates(iReplicates);
        // ### CAUTION ### : need to implement the exact same output as DAISIE_sim
        // how to combine the multiple data types? and which types btw?

        // set max species ID to amount of mainland species
        Archi::setMaxID(iMainSp_n);

        // loop through replicates
        for (int i = 1; i <= iReplicates; ++i) {

            // initialise intermediate archipelago data frame
            Archipelago aAggregArchi(iNumIslands, iAK);

            // run simulation for each mainland sp. seperately --> clade-specific carrying capacity
            for (int j = 0; j < iMainSp_n; ++j) {

                aAggregArchi.addArchi(ArchiDAISIE_core(dAge, 1, vIniPars, iAK, iNumIslands, prng));
                // ### CAUTION: ### :  have to implement that mainland species have different IDs (ID = j)
                    // or: make maxID static so that it increments for all Archipelagos
                    // (shared variable among all objects of class)
            }
            vFinalIslandReplicates[i] = aAggregArchi.aggregateArchi();
        }
        return vFinalIslandReplicates;
    }
    catch (exception &error) {
        std::cerr << "Error: " << error.what();
        exit(1);
    }
}


int main() {

    cout << "hi" << endl;
    vector<double> vPars( {0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12, 50} );
    ArchiDAISIE(5, 50, vPars, 2, 100);
    mt19937_64 prng;
    vector<double> vIni = vPars;
    vIni.pop_back();
    Archi::setMaxID(50);
    Archipelago arch = ArchiDAISIE_core(2, 50, vIni, static_cast<int>(vPars[8]), 3, prng);
    arch.printArchi();
    vector<Isl> test = arch.returnArchi();
    cout << test.empty() << endl;

    return 0;
}//testt
//testt
