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

Archipelago ArchiDAISIE_core(
    const double island_age,
    const int n_mainland_species,
    const vector<double> &vIniPars,
    const int archipelago_carrying_capacity,
    const int n_islands,
    mt19937_64 &prng
)
{
    try {
        // initialise Archipelago data frame and set time and max species ID to initial values
        Archipelago aArchi(n_islands, archipelago_carrying_capacity);
        double dTime = island_age;

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
            vector<double> pLAndGRates = aArchi.calculateAllRates(vIniPars, (int) n_mainland_species, n_islands);

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
            vector<int> vHappening = aArchi.sampleNextEvent(pLAndGRates, prng, (int) n_mainland_species);

            // update the phylogeny
            aArchi.updateArchi(vHappening, vIniPars[1], prng, dTime);
        }
        return aArchi;
    }
    catch (string &str) {
        std::cerr << "Warning: " << str;
        assert(!"should never get here!");
        return Archipelago();
    }
}

vector<vector<Species> > ArchiDAISIE(const double &dAge, const int iMainSp_n, vector<double> &vIniPars,
                                        const unsigned int &iNumIslands, const int iReplicates)
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
        Archipelago::setMaxID(iMainSp_n);

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


void test_island()
{
    {
        const int k{12};
        const Island island(k);
        assert(k ==island.getCarryingCap());
    }
    {
      const int
      Island island(10);
      assert(island.getNspecies() == 0);
      island.pushbackSp(Species(0, -1, 0));
      assert(island.getNspecies() == 1);
    }
}

int main() {

    test_island();
    vector<double> vPars( {0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12, 50} );
    ArchiDAISIE(5, 50, vPars, 2, 100);
    mt19937_64 prng;
    vector<double> vIni = vPars;
    vIni.pop_back();
    Archipelago::setMaxID(50);
    Archipelago arch = ArchiDAISIE_core(2, 50, vIni, static_cast<int>(vPars[8]), 3, prng);

    return 0;
}
