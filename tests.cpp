//
// Created by Sebastian Mader on 26.09.2018.
//

#include "tests.h"

using namespace std;

void test_island()
{
    #ifdef ISLAND_HAS_DEFAULT_CONSTRUCTOR
    {   // A default constructed island has a carrying capacity of zero
        const Island island;
        assert(island.getCarryingCap() == 0);
    }
    #endif
    {   // Setting and getting carrying capacity matches
        // RJCB: a test should have a properly documented purpose;
        // 'checking island' provides no information about the test's intent
        const int k{ 12 };
        const Island island(k);
        assert(k==island.getCarryingCap());
    }
    {   // Adding a species increased the number of species
        // RJCB: a test should only test one thing
        Island island(10);
        assert(island.getNSpecies()==0);
        island.addSpecies(Species());
        assert(island.getNSpecies()==1);
    }
    {   // An added species is an immigrant
        // RJBC: I feel this is weird: there is a function called 'immigrate',
        // so I do not see the use of 'addSpecies'.
        Island island(10);
        island.addSpecies(Species());
        assert(island.getSpecies()[0].readStat() == 'I');
    }
    {   // The immigration of an absent species increases the number of species
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 3.14); //RJCB: insteaf of the int 42, I would prefer a 'Species' or 'SpeciesID' instead
        assert(island.getNSpecies()==1);
        //assert(island.findSpecies(42).readBirth()==3.14);
        //assert(island.findSpecies(42).readStat() == 'I');
    }
    {   // The immigration of an absent species is recognized as an immigrant
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.findSpecies(SpeciesID(42)).readStat() == 'I');
    }
    {   // An immigrant gets its immigration time stored
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 3.14); //RJCB: insteaf of the int 42, I would prefer a 'Species' or 'SpeciesID' instead
        assert(island.findSpecies(SpeciesID(42)).readBirth()==3.14);
    }
    {   // When a same species immigrates twice, the last imigration time is stored
        // RJCB: I would enjoy a class called 'Time' (with only a simple double),
        // so it can be documented if the time here is either
        //  (1) after the start of the process
        //  (2) before the present
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 6.28);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.findSpecies(SpeciesID(42)).readBirth()==3.14); //RJCB: Unsure, see above
        assert(island.getNSpecies()==1);
    }
    {   // Extinction decreases the number of species
        // RJCB: this is the simplest test possible to test extinction
        Island island(10);
        island.immigrate(SpeciesID(42), 6.28);
        assert(island.getNSpecies()==1);
        island.goExtinct(SpeciesID(42));
        assert(island.getNSpecies()==0);
    }
    {   // Extinction decreases the number of species
        // RJCB: this is not the simplest test possible
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 6.28);
        island.immigrate(SpeciesID(1), 3.14);
        assert(island.getNSpecies()==2);
        island.goExtinct(SpeciesID(42));
        assert(island.getNSpecies()==1);
        // RJCB: implementation detail (-1 is not a number inspired by the real world).
        // Use
        //
        //   assert(island.hasSpecies(42));
        //
        // instead. I would expect 'bool Island::hasSpecies(const SpeciesID&) const noexcept' to be present.
        //
        assert(island.findPos(SpeciesID(42)) == -1); //RJCB: this is a seperate test, see the one below
    }
    #ifdef SPECIES_CANNOT_BE_FOUND_ON_EMPTY_ISLAND
    {   // Species cannot be found on empty island
        Island island(10);
        assert(!island.hasSpecies(SpeciesID(42)));
    }
    {   // Species cannot be found on island after immigration
        Island island(10);
        assert(!island.hasSpecies(SpeciesID(42)));
        island.immigrate(SpeciesID(42), 6.28);
        assert(island.hasSpecies(SpeciesID(42)));
    }
    #endif //SPECIES_CANNOT_BE_FOUND_ON_EMPTY_ISLAND
    #ifdef EXTINCTION_OF_ABSENT_SPECIES_THROWS_AN_EXCEPTION
    // Extinction of absent species throws an exception
    {
        Island island(1);
        try
        {
            island.goExtinct(SpeciesID(42));
            assert(!"Should not get here"); //!OCLINT accepted idiom
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "absent species cannot go extinct");
        }
    }
    #endif // EXTINCTION_OF_ABSENT_SPECIES_THROWS_AN_EXCEPTION
    #ifdef IMMIGRATION_OVERFLOODING_CARRYING_CAPACITY_MUST_THROW
    // Species cannot immigrate if island which has as much species as its carrying capacity
    {
        Island island(0);
        try
        {
            island.immigrate(SpeciesID(42), 3.14);
            assert(!"Should not get here"); //!OCLINT accepted idiom
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "cannot immigrate to island with a saturated number of species");
        }
    }
    #endif // IMMIGRATION_OVERFLOODING_CARRYING_CAPACITY_MUST_THROW
    {   // testing speciesID + speciation
        //RJCB: test one thing at a time, seperately. The test discription
        //is already vague, start by writing a concrete description of
        //what you are testing, like I did above
        Island island(10);
        const int n_mainlandSpecies = 50;
        SpeciesID maxSpeciesID(n_mainlandSpecies);
        island.immigrate(SpeciesID(1), 3.14);
        island.immigrate(SpeciesID(42), 3.01);
        island.goExtinct(SpeciesID(42));
        island.speciateAna(SpeciesID(1), maxSpeciesID);
        assert(island.findSpecies(maxSpeciesID).readStat() == 'A');
        assert(maxSpeciesID.getMaxSpeciesID()
            == n_mainlandSpecies+1);
        island.immigrate(SpeciesID(42), 2.56);
        island.speciateClado(SpeciesID(42), 2.50, maxSpeciesID);
        assert(island.findSpecies(maxSpeciesID).readStat() == 'C');
        assert(maxSpeciesID.getMaxSpeciesID()
                == n_mainlandSpecies+3);
        assert(island.getNSpecies() == 3);
    }
    {   // testing calculating of rates
        Island island1(10);
        Island island2(20);
        double sumLog = getLogGrowth(island1) + getLogGrowth(island2);
        const int n_mainlandSpecies = 50;
        const int n_islands = 2;
        vector<double> islPars = { 0.1, 0.5, 0.2, 0.2, 0.15 };
        island1.calculateIslRates(islPars, n_mainlandSpecies,
                n_islands, sumLog);
        double sumRates1 = extractSumOfRates(island1);
        const int n_alive1 = island1.getNSpecies();
        const int islandK1 = island1.getCarryingCap();
        assert(1-static_cast<double>(n_alive1)/islandK1 == getLogGrowth(island1));
        const double immiRate1 = max(0.0, islPars[0] * n_mainlandSpecies
                * getLogGrowth(island1) / n_islands);
        assert(sumRates1 == immiRate1);
    }
    {   // testing sampling of local event
        Island island1(10);
        Island island2(20);
        double sumLog = getLogGrowth(island1) +getLogGrowth(island2);
        const int n_mainlandSpecies = 50;
        const int n_islands = 2;
        vector<double> islPars = { 0.1, 0.5, 0.2, 0.2, 0.15 };
        island1.calculateIslRates(islPars, n_mainlandSpecies,
                n_islands, sumLog);
        double sumRates1 = extractSumOfRates(island1);
        const double immiRate1 = max(0.0, islPars[0] * n_mainlandSpecies
                * getLogGrowth(island1) / n_islands);
        assert(sumRates1 == immiRate1);
        mt19937_64 prng;
        // event_type event = island1.sampleLocalEvent(prng);
    }
}

void test_archi()
{
    {  // test constructor and the correct implementation of carrying capacity
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        Archipelago archi(n_islands, archi_carryingCap);
        assert(archi.getNSpecies()==0);
        assert(archi.getCarryingCap()==archi_carryingCap);
        assert(archi.getIslands()[0].getCarryingCap()==archi_carryingCap/n_islands);
        assert(archi.getIslands()[1].getNSpecies()==0);
    }
    {  // testing the calculation of event rates and initialisation of rate vectors
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        Archipelago archi(n_islands, archi_carryingCap);
        assert(archi.getGlobalRates().empty());
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        archi.calculateAllRates(pars, n_mainland, n_islands);
        assert(archi.getGlobalRates().size()==3);
        vector<Island> archiCopy = archi.getIslands();
        assert(extractSumOfRates(archiCopy[0])>0);
        assert(extractSumOfRates(archiCopy[1])>0);
        assert(archiCopy[0].getLocalRates().size()==5);
        assert(archiCopy[1].getLocalRates().size()==5);
    }
    {  // testing the sampling of next event
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        assert(archi.getGlobalRates().empty());
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        archi.calculateAllRates(pars, n_mainland, n_islands);
        event_type event = archi.sampleNextEvent(prng);
        assert(getEventNum(event) >= 0);
    }
}
