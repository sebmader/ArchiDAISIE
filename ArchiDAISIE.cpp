//
// Created by Sebastian Mader on 07.11.2018.
//

#include "ArchiDAISIE.h"

using namespace std;
namespace fs = experimental::filesystem;

Archipelago ArchiDAISIE_core(const double& islandAge,
        const vector<SpeciesID>& mainSpeciesIDs,
        const vector<double>& initialParameters,
        const int islCarryingCap,
        const int n_islands,
        mt19937_64& prng,
        SpeciesID& maxSpeciesID,
        STTtable& STT)
{
    try {
        // initialise Archipelago data frame and
        // set time to island age (= emergence time of island)
        Archipelago archi(n_islands, islCarryingCap);
        double timeNow = islandAge;
        const int n_mainland = static_cast<int>(mainSpeciesIDs.size());

        // start looping through time
        for (;;) {

            // calculate the rates of events
            archi.calculateAllRates(initialParameters, n_mainland, n_islands);

            // draw time interval to next event
            const vector<double> globalRates = archi.getGlobalRates();
            double sumOfRates = globalRates[0] + globalRates[1] + globalRates[2];
            for (auto& island : archi.getIslands()) {
                sumOfRates += extractSumOfRates(island);
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
            event_type nextEvent = archi.sampleNextEvent(prng);

            // update the phylogeny
            archi.doNextEvent(nextEvent, initialParameters[1], prng, timeNow,
                    maxSpeciesID, mainSpeciesIDs);

            // update STT
            STT.updateFullSTTtable(archi,timeNow);
        }
        return archi;
    }
    catch (string &str) {
        cerr << "Warning: " << str;
        assert(!"should never get here!");  //!OCLINT
        return Archipelago();
    }
}

vector<Island> ArchiDAISIE(const double& islandAge,
        const int n_mainlandSpecies,
        vector<double> initialParameters,
        const int n_islands,
        const int replicates,
        const string& output_dir,
        const int n_timeSlicesSTT)
{
    try {
        // check given parameters
        if (islandAge < 0)
            throw logic_error("Age has to be higher than zero.");
        if (n_mainlandSpecies < 0)
            throw logic_error("Simulation needs at least one mainland species.");
        if(n_islands < 0)
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
        const int islCarryingCap = static_cast<int>(initialParameters[8]);
        initialParameters.pop_back();

        // initialise main data frame
        vector<Island> islandReplicates((unsigned)replicates);
        // ### CAUTION ### : need to implement the exact same output as DAISIE_sim
        // how to combine the multiple data types? and which types btw?

        // output of simulation parameters & seed
        fs::create_directories(output_dir);
        ofstream ofsSimData(output_dir + "/sim_data.txt");
        if (!ofsSimData.is_open())
            throw runtime_error("Unable to open output file.");
        ofsSimData << "Seed: " << seed << '\n'
                   << "Island age: " << islandAge << '\n'
                   << "N Islands: " << n_islands << '\n'
                   << "N Mainland: " << n_mainlandSpecies << '\n'
                   << "Parameters: " << '\n';
        ofsSimData << "immigration: " << initialParameters[0] << '\n'
                   << "migration: " << initialParameters[1] << '\n'
                   << "cladogenesis_local: " << initialParameters[2] << '\n'
                   << "anagenesis_local: " << initialParameters[3] << '\n'
                   << "extinction_local: " << initialParameters[4] << '\n'
                   << "cladogenesis_global: " << initialParameters[5] << '\n'
                   << "anagenesis_global: " << initialParameters[6] << '\n'
                   << "extinction_global: " << initialParameters[7] << '\n'
                   << "island K: " << initialParameters[8] << '\n';
        ofsSimData.close();

        // loop through replicates
        for (int rep = 0; rep < replicates; ++rep) {

            cout << "Simulating replicate " << rep+1 << '\n';

            // initialise intermediate STT data frame
            vector<STTtable> sttPerColoniser((unsigned)n_mainlandSpecies,
                    STTtable(1,STT(islandAge)));
            for (auto& sttTable : sttPerColoniser) {
                assert(sttTable.getSTTtable()[0].getTime()==islandAge);
            }

            // initialise intermediate archipelago data frame
            Archipelago fullArchi(n_islands, islCarryingCap);

            // initialise max species ID
            SpeciesID maxSpeciesID(n_mainlandSpecies);

            // run simulation for each mainland sp. separately -> clade-specific carrying capacity
            for (int mainSp = 0; mainSp < n_mainlandSpecies; ++mainSp) {

                const vector<SpeciesID> mainlandSpecies(1,SpeciesID(mainSp));
                fullArchi.addArchi(ArchiDAISIE_core(islandAge, mainlandSpecies,
                        initialParameters, islCarryingCap,
                        n_islands, prng, maxSpeciesID,
                        sttPerColoniser[mainSp]));
            }
            islandReplicates[rep] = fullArchi.makeArchiTo1Island();

            const STTtable fullSTT = mergeSTTtables(sttPerColoniser, n_timeSlicesSTT);
            // output of merged STT per replicate to file
            ofstream ofsSTT(output_dir + "/rep_" + to_string(rep+1) + "_STT.txt");
            if (!ofsSTT.is_open())
                throw runtime_error("Unable to open output file.");
            ofsSTT << fullSTT;
            ofsSTT.close();
            // output of branching data table per replicate
            ofstream ofsBranching(output_dir + "/rep_" + to_string(rep+1) + "_branching.txt");
            if (!ofsBranching.is_open())
                throw runtime_error("Unable to open output file.");
            outputBranching(islandReplicates[rep],ofsBranching);
            ofsBranching.close();
        }
        return islandReplicates;
    }
    catch (exception &error) {
        cerr << "Error: " << error.what() << '\n';
        exit(1);
    }
}