//
// Created by Sebastian Mader on 26.09.2018.
//

#include "test_island.h"

using namespace std;

void test_island()
{
    {   // testing basics
        const int k{ 12 };
        const Island island(k);
        assert(k==island.getCarryingCap());
    }
    {
        Island island(10);
        assert(island.getNAllSpecies()==0);
        island.addSpecies(Species(0, 0, 0));
        assert(island.getNAllSpecies()==1);
        assert(island.getNSpeciesAlive()==1);
    }
    {
        Island island(10);
        assert(island.getNAllSpecies()==0);
        island.immigrate(42, 3.14);
        assert(island.getNAllSpecies()==1);
        assert(island.findSpecies(42).readBirth()==3.14);
        assert(island.getNSpeciesAlive()==1);
    }
    {
        Island island(10);
        assert(island.getNAllSpecies()==0);
        island.immigrate(42, 6.28);
        island.immigrate(42, 3.14);
        assert(island.findSpecies(42).readBirth()==3.14);
        assert(island.getNAllSpecies()==1);
        assert(island.getNSpeciesAlive()==1);
    }
    {
        Island island(10);
        assert(island.getNAllSpecies()==0);
        assert(island.getNSpeciesAlive()==0);
        island.immigrate(42, 6.28);
        island.immigrate(1, 3.14);
        assert(island.getNAllSpecies()==2);
        assert(island.getNSpeciesAlive()==2);
        island.goExtinct(42, 7.30);
        assert(island.getNAllSpecies()==1);
        assert(island.findPos(42) == -1);
        assert(island.getNSpeciesAlive()==1);
    }
    {
        Island island(10);
        const int n_mainlandSpecies = 50;
        SpeciesID maxSpeciesID(n_mainlandSpecies);
        island.immigrate(1, 3.14);
        island.immigrate(42, 3.01);
        island.goExtinct(42, 2.95);
        island.speciateAna(1, 2.80, maxSpeciesID);
        assert(maxSpeciesID.getMaxSpeciesID()
            == n_mainlandSpecies+1);
        island.immigrate(42, 2.56);
        island.speciateClado(42, 2.50, maxSpeciesID);
        assert(maxSpeciesID.getMaxSpeciesID()
                == n_mainlandSpecies+3);
        assert(island.getNAllSpecies() == 3);
        assert(island.getNSpeciesAlive() == 3);
    }
    {
        Island island1(10);
        Island island2(20);
        double sumLog = island1.returnLogGrowth() + island2.returnLogGrowth();
        const int n_mainlandSpecies = 50;
        const int n_islands = 2;
        vector<double> islPars = { 0.1, 0.5, 0.2, 0.2, 0.15 };
        island1.calculateIslRates(islPars, n_mainlandSpecies,
                n_islands, sumLog);
        double sumRates1 = island1.extractSumOfRates();
        const int n_alive1 = island1.getNSpeciesAlive();
        const int islandK1 = island1.getCarryingCap();
        assert(1-static_cast<double>(n_alive1)/islandK1 == island1.returnLogGrowth());
        const double immiRate1 = max(0.0, islPars[0] * n_mainlandSpecies
                * island1.returnLogGrowth() / n_islands);
        assert(sumRates1 == immiRate1);
    }
    {
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
        mersenne_twister_engine<uint_fast64_t, 64, 312, 156, 31, 0xb5026f5aa96619e9ULL, 29, 0x5555555555555555ULL, 17, 0x71d67fffeda60000ULL, 37, 0xfff7eee000000000ULL, 43, 6364136223846793005ULL> prng;
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
        mersenne_twister_engine<uint_fast64_t, 64, 312, 156, 31, 0xb5026f5aa96619e9ULL, 29, 0x5555555555555555ULL, 17, 0x71d67fffeda60000ULL, 37, 0xfff7eee000000000ULL, 43, 6364136223846793005ULL> prng;
        vector<int> happening = island1.sampleLocalEvent(prng, n_mainlandSpecies);
        island1.immigrate(happening[1], 3.8);
        happening = island1.sampleLocalEvent(prng, n_mainlandSpecies);
        vector<double> logGrowthTerms = { island1.returnLogGrowth(), island2.returnLogGrowth() };
        const int destinationIsl = island1.drawMigDestinationIsland(0,
                                            logGrowthTerms, islPars[1], prng);
        assert(destinationIsl == 1);
        island2.migrate(island1.findSpecies(happening[1]));
        island1.speciateClado(happening[1],3.5, maxSpeciesID);
        sumLogWO1 = island2.returnLogGrowth();
        island1.calculateIslRates(islPars, n_mainlandSpecies, n_islands, sumLogWO1);
        happening = island1.sampleLocalEvent(prng, n_mainlandSpecies);
        island2.immigrate(50, 3.2);
        double sumLogWO2 = island1.returnLogGrowth();
        island2.calculateIslRates(islPars, n_mainlandSpecies, n_islands, sumLogWO2);
        island1.immigrate(50, 2.8);
        island1.printIsland();
        island2.printIsland();
        island1.consolidateIslands(island2);
        island1.printIsland();
    }
}