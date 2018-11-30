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
#include <filesystem>


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

        test_other_functions();
        test_speciesID();
        test_species();
        test_island();
        test_archi();
        test_STT();
        test_STTtable();

        string inputFileName = "test.txt";

        const int count = countCSVFileRows(inputFileName);
        int it = 1;
        while(it < count) {
            ++it;
            ifstream ifs(inputFileName);
            vector<string> parameters = readParameterRowCSV(ifs,it);
            assert(parameters.size() == 14);
            const string simName = parameters[0];
            const double immi = stod(parameters[1]);
            const double mig = stod(parameters[2]);
            const double clado_l = stod(parameters[3]);
            const double ana_l = stod(parameters[4]);
            const double ext_l = stod(parameters[5]);
            const double clado_g = stod(parameters[6]);
            const double ana_g = stod(parameters[7]);
            const double ext_g = stod(parameters[8]);
            const int kPerIsl = stoi(parameters[9]);
            const int n_isl = stoi(parameters[10]);
            const double archi_age = stod(parameters[11]);
            const int n_mainlandSp = stoi(parameters[12]);
            const int replicates = stoi(parameters[13]);
            vector<double> initialPars = { immi,mig,clado_l,ana_l,ext_l,clado_g,
                                           ana_g,ext_g };
            string output_dir(fs::current_path().fs::path::parent_path().string()
                    + "/sims/" + simName);
            cout << "Simulation " << it-1 << " / " << count-1 << '\n';
            ArchiDAISIE(archi_age,
                    initialPars,
                    n_mainlandSp,
                    kPerIsl,
                    n_isl,
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
