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
#include <exception>
#include <experimental/filesystem>


using namespace std;
namespace fs = experimental::filesystem;

// ------------ CLASS/FUNCTION DEFINITIONS ------------ //

#include "tests.h"
#include "free_functions.h"
#include "SpeciesID.h"
#include "Species.h"
#include "Island.h"
#include "Archipelago.h"
#include "ArchiDAISIE.h"

// ------------ ArchiDAISIE FUNCTIONS ------------ //

int main() {

    try {
        /*
        test_other_functions();
        test_speciesID();
        test_species();
        test_island();
        test_archi();
        test_STT();
        test_STTtable();
        */
        const int n_sims = 5;
        const int islandAge = 10;
        const int n_mainlandSp = 100;
        const int n_islands = 2;
        const int replicates = 10;
        const int islCarryingCap = 10;

        for (int i = 1; i <= n_sims; ++i) {
            string output_dir(fs::current_path().fs::path::parent_path().string()
                        + "/test_sims/sim_" + to_string(i));
            const vector<double> vPars( {0.1, 0.3, 0.2, 0.12, 0.2, 0.2, 0.1, 0.12,
                                         islCarryingCap} );
            ArchiDAISIE(islandAge,
                    n_mainlandSp,
                    vPars,
                    n_islands,
                    replicates,
                    output_dir,
                    25);
        }
    }
    catch (exception &error) {
        cerr << "Main_Error: " << error.what() << '\n';
        exit(1);
    }

    return 0;
}
