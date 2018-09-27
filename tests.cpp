//
// Created by Sebastian Mader on 26.09.2018.
//

#include "tests.h"

using namespace std;

void test_island()
{
    {   // testing constructor
        const int k{ 12 };
        const Island island(k);
        assert(k==island.getCarryingCap());
    }
    {   // testing adding species
        Island island(10);
        assert(island.getNSpecies()==0);
        island.addSpecies(Species(0, 0, 0));
        assert(island.getNSpecies()==1);
        assert(island.returnIsland()[0].readStat() == 'I');
    }
    {   // testing immigration
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(42, 3.14);
        assert(island.getNSpecies()==1);
        assert(island.findSpecies(42).readBirth()==3.14);
        assert(island.findSpecies(42).readStat() == 'I');
    }
    {   // testing re-immigration
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(42, 6.28);
        island.immigrate(42, 3.14);
        assert(island.findSpecies(42).readBirth()==3.14);
        assert(island.getNSpecies()==1);
    }
    {   // testing extinction
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(42, 6.28);
        island.immigrate(1, 3.14);
        assert(island.getNSpecies()==2);
        island.goExtinct(42);
        assert(island.getNSpecies()==1);
        assert(island.findPos(42) == -1);
    }
    {   // testing maxSpeciesID + speciation
        Island island(10);
        const int n_mainlandSpecies = 50;
        SpeciesID maxSpeciesID(n_mainlandSpecies);
        island.immigrate(1, 3.14);
        island.immigrate(42, 3.01);
        island.goExtinct(42);
        island.speciateAna(1, maxSpeciesID);
        assert(island.findSpecies(n_mainlandSpecies+1).readStat() == 'A');
        assert(maxSpeciesID.getMaxSpeciesID()
            == n_mainlandSpecies+1);
        island.immigrate(42, 2.56);
        island.speciateClado(42, 2.50, maxSpeciesID);
        assert(island.findSpecies(n_mainlandSpecies+2).readStat() == 'C');
        assert(island.findSpecies(n_mainlandSpecies+3).readStat() == 'C');
        assert(maxSpeciesID.getMaxSpeciesID()
                == n_mainlandSpecies+3);
        assert(island.getNSpecies() == 3);
    }
    {   // testing calculating of rates
        Island island1(10);
        Island island2(20);
        double sumLog = island1.returnLogGrowth() + island2.returnLogGrowth();
        const int n_mainlandSpecies = 50;
        const int n_islands = 2;
        vector<double> islPars = { 0.1, 0.5, 0.2, 0.2, 0.15 };
        island1.calculateIslRates(islPars, n_mainlandSpecies,
                n_islands, sumLog);
        double sumRates1 = island1.extractSumOfRates();
        const int n_alive1 = island1.getNSpecies();
        const int islandK1 = island1.getCarryingCap();
        assert(1-static_cast<double>(n_alive1)/islandK1 == island1.returnLogGrowth());
        const double immiRate1 = max(0.0, islPars[0] * n_mainlandSpecies
                * island1.returnLogGrowth() / n_islands);
        assert(sumRates1 == immiRate1);
    }
    {   // testing sampling of local event
        Island island1(10);
        Island island2(20);
        double sumLog = island1.returnLogGrowth() + island2.returnLogGrowth();
        const int n_mainlandSpecies = 50;
        const int n_islands = 2;
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
    }
    {   // testing all
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
        happening = island1.sampleLocalEvent(prng, n_mainlandSpecies);
        vector<double> logGrowthTerms = { island1.returnLogGrowth(), island2.returnLogGrowth() };
        const int destinationIsl = island1.drawMigDestinationIsland(0,
                                            logGrowthTerms, islPars[1], prng);
        assert(destinationIsl == 1);
        island2.migrate(island1.findSpecies(happening[1]), 3.6);
        island1.speciateClado(happening[1],3.5, maxSpeciesID);
        sumLogWO1 = island2.returnLogGrowth();
        island1.calculateIslRates(islPars, n_mainlandSpecies, n_islands, sumLogWO1);
        happening = island1.sampleLocalEvent(prng, n_mainlandSpecies);
        island2.immigrate(50, 3.2);
        double sumLogWO2 = island1.returnLogGrowth();
        island2.calculateIslRates(islPars, n_mainlandSpecies, n_islands, sumLogWO2);
        island1.immigrate(50, 2.8);
        island2.migrate(island1.findSpecies(50), 2.5);
        island2.immigrate(23, 2.2);
        island2.speciateClado(23, 2.0, maxSpeciesID);
        island1.migrate(island2.findSpecies(maxSpeciesID.getMaxSpeciesID()),1.74);
        island1.printIsland();
        island2.printIsland();
        island1.consolidateIslands(island2);
        island1.printIsland();
    }
}