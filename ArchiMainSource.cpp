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

using namespace std;

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

    test_other_functions();
    test_STT();
    test_STTtable();
    test_speciesID();
    test_species();
    test_island();
    test_archi();

    const int islandAge = 4;
    const int n_mainlandSp = 100;
    const int n_islands = 2;
    const int replicates = 1;
    const vector<double> vPars( {0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12, 50} );
    const vector<Island> archipelago = ArchiDAISIE(islandAge,n_mainlandSp,vPars,n_islands,replicates);

    int i = 0;
    for (auto& z : archipelago) {
        ++i;
        cout << "Island " << i << '\n';
        cout << "BirT" << '\t' << "ParID" << '\t'
             << "SpID" << '\t' << "Status" << '\t'
             << "Mig" << '\t'  << "AncesT" << '\t'
             << "ColoT" << '\t' << "CladoStac" << '\n';
        z.printIsland();
    }

    return 0;
}
