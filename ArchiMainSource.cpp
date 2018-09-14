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

// ------------ FUNCTION DEFINITIONS ------------ //


// ------------ CLASS DEFINITIONS ------------ //

#include "Species.h"
#include "Island.h"
#include "Archipelago.h"

// ------------ ArchiDAISIE FUNCTIONS ------------ //

Archipelago ArchiDAISIE_core(const double islandAge,
        const int n_mainlandSpecies,
        const vector<double>& initialParameters,
        const int archiCarryingCap,
        const int n_islands,
        mt19937_64& prng,
        SpeciesID& maxSpeciesID)
{
    try {
        // initialise Archipelago data frame and
        // set time to island age (= emergence time of island)
        Archipelago aArchi(n_islands, archiCarryingCap);
        double timeNow = islandAge;

        // start looping through time
        for (;;) {

            // calculate the rates of events
            vector<double> localGlobalRates = aArchi.calculateAllRates(
                    initialParameters, n_mainlandSpecies, n_islands);
            assert(localGlobalRates.size() == 2);

            // draw time interval to next event
            const double sumOfRates = localGlobalRates[0] + localGlobalRates[1];
            if (sumOfRates <= 0)
                throw runtime_error("Event rate is zero or below. "
                                    "No event can be drawn.\n");

            // draw time to next event and take it off the actual time
            exponential_distribution<double> waitingTime(sumOfRates);
            const double dT = waitingTime(prng);
            timeNow -= dT;
            if (timeNow <= 0)  // when timeNow passes the present, stop simulation
                break;

            // sample which event happens
            vector<int> vHappening = aArchi.sampleNextEvent(localGlobalRates,
                    prng, n_mainlandSpecies);

            // update the phylogeny
            aArchi.doNextEvent(vHappening, initialParameters[1], prng, timeNow,
                    maxSpeciesID);
        }
        return aArchi;
    }
    catch (string &str) {
        std::cerr << "Warning: " << str;
        assert(!"should never get here!");
        return Archipelago();
    }
}

vector<vector<Species> > ArchiDAISIE(const double &islandAge,
        const int n_mainlandSpecies,
        vector<double> &initialParameters,
        const unsigned int &n_islands,
        const int replicates)
{   // ### CAUTION ### : output unclear !!
    try {
        // check given parameters
        if (islandAge <= 0)
            throw logic_error("Age has to be higher than zero.");
        if (n_mainlandSpecies <= 0)
            throw logic_error("Simulation needs at least one mainland species.");
        if(n_islands <= 0)
            throw logic_error("Simulation needs at least one island.");
        if(initialParameters.size() != 9)
            throw logic_error("Provide 9 parameter values.");
        if(initialParameters[0] <= 0)
            throw logic_error("Rate of colonisation is zero or below."
                              " The island cannot be colonised.");

        // declare and seed PRNG with system clock
        mt19937_64 prng;

        chrono::high_resolution_clock::time_point tp =
                chrono::high_resolution_clock::now();
        const unsigned seed = static_cast<unsigned>(
                tp.time_since_epoch().count());
        // ### CAUTION: ###   should print/output seed here
        prng.seed(seed);

        // order of parameters (input):
        // gam_i, gam_m, lamb_cl, lamb_al, mu_l,
        // lamb_cg, lamb_ag, mu_g, archi-wide K
        const int archiCarryingCap = static_cast<int>(initialParameters[8]);
        initialParameters.pop_back();

        // initialise main data frame
        vector< vector<Species> > vFinalIslandReplicates(replicates);
        // ### CAUTION ### : need to implement the exact same output as DAISIE_sim
        // how to combine the multiple data types? and which types btw?

        // loop through replicates
        for (int i = 1; i <= replicates; ++i) {

            // initialise intermediate archipelago data frame
            Archipelago aAggregArchi(n_islands, archiCarryingCap);

            // initialise max species ID
            SpeciesID maxSpeciesID(n_mainlandSpecies);

            // run simulation for each mainland sp. separately -> clade-specific carrying capacity
            for (int j = 0; j < n_mainlandSpecies; ++j) {

                aAggregArchi.addArchi(ArchiDAISIE_core(islandAge, 1, initialParameters,
                        archiCarryingCap, n_islands, prng, maxSpeciesID));
            }
            vFinalIslandReplicates[i] = aAggregArchi.makeArchiTo1Island();
        }
        return vFinalIslandReplicates;
    }
    catch (exception &error) {
        std::cerr << "Error: " << error.what();
        exit(1);
    }
}


void test_island()
{
    {
        const int k{ 12 };
        const Island island(k);
        assert(k==island.getCarryingCap());
    }
    {
        Island island(10);
        assert(island.getNAllSpecies()==0);
        island.addSpecies(Species(0, 0, 0));
        assert(island.getNAllSpecies()==1);
        assert(island.getNSpeciesAlive()==1);
    }
    {
        Island island(10);
        assert(island.getNAllSpecies()==0);
        island.immigrate(42, 3.14);
        assert(island.getNAllSpecies()==1);
        assert(island.findSpecies(42).readBirth()==3.14);
        assert(island.getNSpeciesAlive()==1);
    }
    {
        Island island(10);
        assert(island.getNAllSpecies()==0);
        island.immigrate(42, 3.14);
        island.immigrate(42, 6.28);
        assert(island.getNAllSpecies()==1);
        assert(island.getNSpeciesAlive()==1);
    }
    {
        Island island(10);
        assert(island.getNAllSpecies()==0);
        island.immigrate(42, 3.14);
        island.immigrate(42, 6.28);
        assert(island.findSpecies(42).readBirth()==6.28);
    }
    {
        Island island(10);
        assert(island.getNAllSpecies()==0);
        assert(island.getNSpeciesAlive()==0);
        island.immigrate(1, 3.14);
        island.immigrate(42, 6.28);
        assert(island.getNAllSpecies()==2);
        assert(island.getNSpeciesAlive()==2);
        island.goExtinct(42, 7.30);
        assert(island.getNAllSpecies()==2);
        assert(!island.findSpecies(42).isExtant());
        assert(island.getNSpeciesAlive()==1);
    }
    {
        Island island(10);
        const int n_mainlandSpecies = 50;
        SpeciesID maxSpeciesID(n_mainlandSpecies);
        island.immigrate(1, 3.14);
        island.immigrate(42, 3.01);
        island.goExtinct(42, 2.95);
        island.speciateAna(1, 2.80, maxSpeciesID);
        assert(maxSpeciesID.getMaxSpeciesID()
            == n_mainlandSpecies+1);
        island.immigrate(42, 2.56);
        island.speciateClado(42, 2.50, maxSpeciesID);
        assert(maxSpeciesID.getMaxSpeciesID()
                == n_mainlandSpecies+3);
        assert(island.getNAllSpecies() == 6);
        assert(island.getNSpeciesAlive() == 3);
    }
    {
        Island island1(10);
        Island island2(20);
        double sumLog = island1.returnLogGrowth() + island2.returnLogGrowth();
        const int n_mainlandSpecies = 50;
        const int n_islands = 2;
        vector<double> islPars = { 0.1, 0.5, 0.2, 0.2, 0.15 };
        island1.calculateIslRates(islPars, n_mainlandSpecies,
                n_islands, sumLog);
        double sumRates1 = island1.extractSumOfRates();
        const int n_alive1 = island1.getNSpeciesAlive();
        const int islandK1 = island1.getCarryingCap();
        assert(1-static_cast<double>(n_alive1)/islandK1 == island1.returnLogGrowth());
        const double immiRate1 = max(0.0, islPars[0] * n_mainlandSpecies
                * island1.returnLogGrowth() / n_islands);
        assert(sumRates1 == immiRate1);
    }
    {
        Island island1(10);
        Island island2(20);
        double sumLog = island1.returnLogGrowth() + island2.returnLogGrowth();
        const int n_mainlandSpecies = 50;
        const int n_islands = 2;
        SpeciesID maxSpeciesID(n_mainlandSpecies);
        vector<double> islPars = { 0.1, 0.5, 0.2, 0.2, 0.15 };
        island1.calculateIslRates(islPars, n_mainlandSpecies,
                n_islands, sumLog);
        double sumRates1 = island1.extractSumOfRates();
        const double immiRate1 = max(0.0, islPars[0] * n_mainlandSpecies
                * island1.returnLogGrowth() / n_islands);
        assert(sumRates1 == immiRate1);
        mt19937_64 prng;
        vector<int> happening = island1.sampleLocalEvent(prng, n_mainlandSpecies);
        assert(happening.size() == 2);
        assert(happening[0] == 0); // has to be immigration
        assert(happening[1] <= n_mainlandSpecies);
        island1.immigrate(happening[1], 3.8);
        Species sp = island1.findSpecies(happening[1]);
        assert(sp.readSpID() == happening[1]);
        assert(island1.findPos(sp.readSpID() == 0));
        assert(island1.returnIsland().size() == 1);
        island1.printIsland();
    }
    {
        Island island1(10);
        Island island2(20);
        double sumLogWO1 = island2.returnLogGrowth();
        const int n_mainlandSpecies = 50;
        const int n_islands = 2;
        SpeciesID maxSpeciesID(n_mainlandSpecies);
        vector<double> islPars = { 0.1, 0.5, 0.2, 0.2, 0.15 };
        island1.calculateIslRates(islPars, n_mainlandSpecies,
                n_islands, sumLogWO1);
        mt19937_64 prng;
        vector<int> happening = island1.sampleLocalEvent(prng, n_mainlandSpecies);
        island1.immigrate(happening[1], 3.8);
        vector<double> logGrowthTerms = { island1.returnLogGrowth(), island2.returnLogGrowth() };
        const int destinationIsl = island1.drawMigDestinationIsland(0,
                                            logGrowthTerms, islPars[1], prng);
        assert(destinationIsl == 1);
    }
}

int main() {

    test_island();
    exit(0);
    vector<double> vPars( {0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12, 50} );
    ArchiDAISIE(5, 50, vPars, 2, 100);
    mt19937_64 prng;
    vector<double> vIni = vPars;
    vIni.pop_back();
    SpeciesID maxSpeciesID(0);
    Archipelago arch = ArchiDAISIE_core(2, 50, vIni, static_cast<int>(vPars[8]),
            3, prng, maxSpeciesID);

    return 0;
}
