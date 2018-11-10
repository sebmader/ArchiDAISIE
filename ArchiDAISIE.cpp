//
// Created by Sebastian Mader on 07.11.2018.
//

#include "ArchiDAISIE.h"

using namespace std;

Archipelago ArchiDAISIE_core(const double& islandAge,
        const int n_mainlandSpecies,
        const std::vector<double>& initialParameters,
        const int archiCarryingCap,
        const int n_islands,
        mt19937_64& prng,
        SpeciesID& maxSpeciesID)
{
    try {
        // initialise Archipelago data frame and
        // set time to island age (= emergence time of island)
        Archipelago archi(n_islands, archiCarryingCap);
        double timeNow = islandAge;

        // start looping through time
        for (;;) {

            // calculate the rates of events
            archi.calculateAllRates(initialParameters, n_mainlandSpecies, n_islands);

            // draw time interval to next event
            const std::vector<double> globalRates = archi.getGlobalRates();
            double sumOfRates = globalRates[0] + globalRates[1] + globalRates[2];
            for (auto& island : archi.getIslands()) {
                sumOfRates += extractSumOfRates(island);
            }
            if (sumOfRates <= 0)
                throw std::runtime_error("Event rate is zero or below. "
                                    "No event can be drawn.\n");

            // draw time to next event and take it off the actual time
            exponential_distribution<double> waitingTime(sumOfRates);
            const double dT = waitingTime(prng);
            timeNow -= dT;
            if (timeNow <= 0)  // when timeNow passes the present, stop simulation
                break;

            // sample which event happens
            event_type nextEvent = archi.sampleNextEvent(prng);

            // update the phylogeny
            archi.doNextEvent(nextEvent, initialParameters[1], prng, timeNow,
                    maxSpeciesID, n_mainlandSpecies);
        }
        return archi;
    }
    catch (std::string &str) {
        std::cerr << "Warning: " << str;
        assert(!"should never get here!");  //!OCLINT
        return Archipelago();
    }
}

std::vector<Island> ArchiDAISIE(const double& islandAge,
        const int n_mainlandSpecies,
        std::vector<double> initialParameters,
        const int n_islands,
        const int replicates)
{   // ### CAUTION ### : output unclear !!
    try {
        // check given parameters
        if (islandAge < 0)
            throw std::logic_error("Age has to be higher than zero.");
        if (n_mainlandSpecies < 0)
            throw std::logic_error("Simulation needs at least one mainland species.");
        if(n_islands < 0)
            throw std::logic_error("Simulation needs at least one island.");
        if(initialParameters.size() != 9)
            throw std::logic_error("Provide 9 parameter values.");
        if(initialParameters[0] <= 0)
            throw std::logic_error("Rate of colonisation is zero or below."
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
        std::vector<Island> islandReplicates(replicates);
        // ### CAUTION ### : need to implement the exact same output as DAISIE_sim
        // how to combine the multiple data types? and which types btw?

        // loop through replicates
        for (int rep = 0; rep < replicates; ++rep) {

            // initialise intermediate archipelago data frame
            Archipelago fullArchi(n_islands, archiCarryingCap);

            // initialise max species ID
            SpeciesID maxSpeciesID(n_mainlandSpecies);

            // run simulation for each mainland sp. separately -> clade-specific carrying capacity
            for (int mainSp = 0; mainSp < n_mainlandSpecies; ++mainSp) {

                fullArchi.addArchi(ArchiDAISIE_core(islandAge, 1, initialParameters,
                        archiCarryingCap, n_islands, prng, maxSpeciesID));
            }
            islandReplicates[rep] = fullArchi.makeArchiTo1Island();
        }
        return islandReplicates;
    }
    catch (std::exception &error) {
        std::cerr << "Error: " << error.what();
        exit(1);
    }
}