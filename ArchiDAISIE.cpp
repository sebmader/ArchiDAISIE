//
// Created by Sebastian Mader on 07.11.2018.
//

#include "ArchiDAISIE.h"

using namespace std;
namespace fs = experimental::filesystem;

Archipelago ArchiDAISIE_core(const double& islandAge, const std::vector<SpeciesID>& mainSpeciesIDs,
        const std::vector<double>& initialParameters, int islCarryingCap, int n_islands, std::mt19937_64& prng,
        SpeciesID& maxSpeciesID, STTtable& STT, int& n_events, int& n_globalEvents)
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

            cout << timeNow << ": ";
            // sample which event happens
            event_type nextEvent = archi.sampleNextEvent(prng);

            ++n_events;
            if(is_global(nextEvent))
                ++n_globalEvents;

            // update the phylogeny
            archi.doNextEvent(nextEvent, initialParameters[1], prng, timeNow,
                    maxSpeciesID, mainSpeciesIDs);

            // update STT
            STT.updateFullSTTtable(archi,timeNow);
        }

        return archi;
    }
    catch (exception &ex) {
        cerr << "ArchiDAISIE_core_error: " << ex.what() << '\n';
        exit(EXIT_FAILURE);
    }
}

vector<Island> ArchiDAISIE(const double& islandAge,
        vector<double> initialParameters,
        int n_mainlandSpecies,
        int kPerIsl,
        int n_islands,
        int replicates,
        const string& output_dir,
        int n_timeSlicesSTT)
{
    try {
        // check given parameters
        if (islandAge < 0)
            throw logic_error("Age has to be higher than zero.");
        if (n_mainlandSpecies < 0)
            throw logic_error("Simulation needs at least one mainland species.");
        if(n_islands < 0)
            throw logic_error("Simulation needs at least one island.");
        if(kPerIsl < 0)
            throw logic_error("The carrying capacity cannot be below 0.");
        if(initialParameters.size() != 8)
            throw logic_error("Provide 9 parameter values.");
        if(initialParameters[0] <= 0)
            throw logic_error("Rate of colonisation is zero or below."
                              " The island cannot be colonised.");

        // order of parameters (input):
        // gam_i, gam_m, lamb_cl, lamb_al, mu_l,
        // lamb_cg, lamb_ag, mu_g

        // declare and seed PRNG with system clock
        mt19937_64 prng;
        chrono::high_resolution_clock::time_point tp =
                chrono::high_resolution_clock::now();
        const unsigned seed = static_cast<unsigned>(
                tp.time_since_epoch().count());
        prng.seed(seed);

        // create output directory
        fs::create_directories(output_dir);

        // initialise main data frame
        vector<Island> islandReplicates((unsigned)replicates);

        double sumNSpecies = 0.0, sumNSpStdDev = 0.0;
        double sumColonisations = 0.0, sumColoStdDev = 0.0;
        double sumGlobalRatio = 0.0, sumGlobalStdDev = 0.0;
        int n_events = 0, n_global_events = 0;

        // loop through replicates
        for (int rep = 0; rep < replicates; ++rep) {

            // cout << "Simulating replicate " << rep+1 << '\n';

            // initialise intermediate STT data frame
            vector<STTtable> sttPerColoniser((unsigned)n_mainlandSpecies,
                    STTtable(1,STT(islandAge)));
            for (auto& sttTable : sttPerColoniser) {
                assert(sttTable.getSTTtable()[0].getTime()==round(islandAge));
            }

            // initialise intermediate archipelago data frame
            Archipelago fullArchi(n_islands, kPerIsl);

            // initialise max species ID
            SpeciesID maxSpeciesID(n_mainlandSpecies);

            // run simulation for each mainland sp. separately -> clade-specific carrying capacity
            for (int mainSp = 1; mainSp <= n_mainlandSpecies; ++mainSp) {

                const vector<SpeciesID> mainlandSpecies(1,SpeciesID(mainSp));
                fullArchi.addArchi(ArchiDAISIE_core(islandAge, mainlandSpecies,
                        initialParameters, kPerIsl,
                        n_islands, prng, maxSpeciesID,
                        sttPerColoniser[mainSp-1], n_events, n_global_events));
            }
            // stats
            sumNSpecies += fullArchi.getNSpecies();
            sumNSpStdDev += fullArchi.getNSpecies() * fullArchi.getNSpecies();
            sumColonisations += fullArchi.getNColonisations();
            sumColoStdDev += fullArchi.getNColonisations() * fullArchi.getNColonisations();
            const double globalRatio = (double)n_global_events/(double)n_events;
            sumGlobalRatio += globalRatio;
            sumGlobalStdDev += globalRatio * globalRatio;

            islandReplicates[rep] = fullArchi.makeArchiTo1Island();
            assert(fullArchi.getNColonisations() == islandReplicates[rep].getNColonisations());

            for (auto& sttTable : sttPerColoniser) {
                assert(sttTable.getSTTtable()[0].getTime()==round(islandAge));
            }
            const STTtable fullSTT = mergeSTTtables(sttPerColoniser, n_timeSlicesSTT);
            // output of merged STT per replicate to file
            ofstream ofsSTT(output_dir + "/rep_" + to_string(rep+1) + "_STT.txt");
            if (!ofsSTT.is_open())
                throw runtime_error("Unable to open STT output file.");
            ofsSTT << fullSTT;
            ofsSTT.close();
            // output of branching data table per replicate
            ofstream ofsBranching(output_dir + "/rep_" + to_string(rep+1) + "_branching.txt");
            if (!ofsBranching.is_open())
                throw runtime_error("Unable to open branching output file.");
            outputBranching(islandReplicates[rep],ofsBranching);
            ofsBranching.close();
        }
        // output of simulation parameters & seed
        ofstream ofsSimData(output_dir + "/sim_data.txt");
        if (!ofsSimData.is_open())
            throw runtime_error("Unable to open simulation data output file.");
        // calculate stats:
        const double meanNSpecies = sumNSpecies/(double)replicates;
        const double nSpeciesStdDev = sqrt((sumNSpStdDev - replicates * meanNSpecies
                * meanNSpecies) / (replicates - 1));
        const double meanColonisations = sumColonisations/(double)replicates;
        const double coloStdDeviation = sqrt((sumColoStdDev - replicates * meanColonisations
                * meanColonisations) / (replicates - 1));
        const double meanGlobalRatio = sumGlobalRatio/(double)replicates;
        const double globalRatioStdDeviation = sqrt((sumGlobalStdDev - replicates * meanGlobalRatio
                * meanGlobalRatio) / (replicates - 1));

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
                   << "island K: " << kPerIsl << '\n';
        ofsSimData << "\nmean number of species: " << setprecision(3)
                   << meanNSpecies << '\n'
                   << "species standard deviation: " << nSpeciesStdDev << '\n';
        ofsSimData << "mean number of colonisations: " << setprecision(3)
                   << meanColonisations << '\n'
                   << "colonisations standard deviation: " << coloStdDeviation << '\n';
        ofsSimData << "mean ratio of global events: " << setprecision(4)
                   << meanGlobalRatio << '\n'
                   << "ratio Global standard deviation: " << globalRatioStdDeviation << '\n';
        ofsSimData.close();

        return islandReplicates;
    }
    catch (exception &error) {
        cerr << "ArchiDAISIE_Error: " << error.what() << '\n';
        exit(EXIT_FAILURE);
    }
}