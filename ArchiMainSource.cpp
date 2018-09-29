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

// ------------ CLASS/FUNCTION DEFINITIONS ------------ //

#include "event_type.h"
#include "Species.h"
#include "Island.h"
#include "Archipelago.h"
#include "tests.h"

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
            aArchi.calculateAllRates(initialParameters, n_mainlandSpecies, n_islands);

            // draw time interval to next event
            vector<double> globalRates = aArchi.getGlobalRates();
            assert(globalRates.size() == 3);
            double sumOfRates = globalRates[0] + globalRates[1] + globalRates[2];
            for (auto& island : aArchi.getIslands()) {
                sumOfRates += island.extractSumOfRates();
            }
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
            event_type nextEvent = aArchi.sampleNextEvent(prng);

            // update the phylogeny
            aArchi.doNextEvent(nextEvent, initialParameters[1], prng, timeNow,
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

int main() {

    test_island();
    test_archi();
    /*
    vector<double> vPars( {0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12, 50} );
    ArchiDAISIE(5, 50, vPars, 2, 100);
    mt19937_64 prng;
    vector<double> vIni = vPars;
    vIni.pop_back();
    SpeciesID speciesID(0);
    Archipelago arch = ArchiDAISIE_core(2, 50, vIni, static_cast<int>(vPars[8]),
            3, prng, speciesID);
    */
    return 0;
}
