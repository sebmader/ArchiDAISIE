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
    }
    {  // test consolidation
        Island island1(10);
        Island island2(20);
        const int n_mainlandSpecies = 50;
        SpeciesID maxSpeciesID(n_mainlandSpecies);
        island1.immigrate(50, 2.8);
        island1.consolidateIslands(island2);
        island2.immigrate(50, 2.7);
        island2.speciateAna(50, maxSpeciesID);
        island2.migrate(island1.findSpecies(50), 2.5);
        island1.migrate(island2.findSpecies(maxSpeciesID.getMaxSpeciesID()), 2.4);
        island2.immigrate(23, 2.2);
        island2.speciateClado(23, 2.0, maxSpeciesID);
        island1.migrate(island2.findSpecies(maxSpeciesID.getMaxSpeciesID()),1.74);
        island1.immigrate(10, 1.9);
        island2.migrate(island1.findSpecies(10), 1.8);
        island1.goExtinct(10);
        island1.migrate(island2.findSpecies(10), 1.5);
        island1.immigrate(1, 1.49);
        island2.immigrate(1, 1.26);
        island1.immigrate(2, 1.3);
        island2.immigrate(2, 1.1);
        island1.immigrate(3, 1.0);
        island2.migrate(island1.findSpecies(3),0.9);
        island2.migrate(island1.findSpecies(3),0.5);
        island1.immigrate(6, 0.4);
        island1.speciateAna(6, maxSpeciesID);
        island2.migrate(island1.findSpecies(maxSpeciesID.getMaxSpeciesID()),0.2);
        island1.consolidateIslands(island2);
    }
}

void test_archi()
{
    {  // test basics
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        Archipelago archi(n_islands, archi_carryingCap);
        assert(archi.getNSpecies() == 0);
        assert(archi.getCarryingCap() == archi_carryingCap);
        assert(archi.returnArchi()[0].getCarryingCap() == archi_carryingCap/n_islands);
        assert(archi.returnArchi()[1].getNSpecies() == 0);
    }
    {
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        Archipelago archi(n_islands, archi_carryingCap);
        assert(archi.getGlobalRates().empty());
        vector<double> pars { 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 } ;
        vector<double> sumLocGlo = archi.calculateAllRates(pars, n_mainland, n_islands);
        assert(archi.getGlobalRates().size() == 3);
        assert(sumLocGlo.size() == 2);
        assert(sumLocGlo[0] == 0);  // global rates == 0 -> no species on archi
        assert(sumLocGlo[1] > 0);  // local rates > 0 -> only immigration can happen
        vector<Island> archiCopy = archi.returnArchi();
        assert(archiCopy[0].extractSumOfRates() > 0);
        assert(archiCopy[1].extractSumOfRates() > 0);
        assert(archiCopy[0].getLocalRates().size() == 5);
        assert(archiCopy[1].getLocalRates().size() == 5);
    }
    {
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        assert(archi.getGlobalRates().empty());
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        vector<double> sumLocGlo = archi.calculateAllRates(pars, n_mainland, n_islands);
        vector<int> happening = archi.sampleNextEvent(sumLocGlo, prng, n_mainland);
        assert(happening.size() == 3); // -> has to be immigration (= local)
        assert(happening[0] == 0 && happening[1] > 0 && happening[2] >= 0);
    }
    {
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        SpeciesID maxSpeciesID(n_mainland);
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        vector<double> sumLocGlo = archi.calculateAllRates(pars, n_mainland, n_islands);
        vector<int> happening = archi.sampleNextEvent(sumLocGlo, prng, n_mainland);
        const double iniMigRate = pars[1];
        archi.doNextEvent(happening, iniMigRate, prng, 3.9, maxSpeciesID);
        vector<int> onWhichIsls = archi.findIsl(happening[1]);
        assert(onWhichIsls.size() == 1); // -> one immigration
        const int inhabitedIslPos = onWhichIsls[0];
        vector<Island> archiCopy = archi.returnArchi();
        const int pos = archiCopy[inhabitedIslPos].findPos(happening[1]);
        assert(pos == 0);
    }
    {
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        SpeciesID maxSpeciesID(n_mainland);
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        const double iniMigRate = 0.1;
        vector<int> happening { 0, 29, 0 };
        archi.doNextEvent(happening, iniMigRate, prng, 4.0, maxSpeciesID);
        vector<int> happening2 { 0, 65, 1 };
        archi.doNextEvent(happening2, iniMigRate, prng, 4.0, maxSpeciesID);
        vector<int> happening3 { 1, 65, 1 };
        archi.doNextEvent(happening3, iniMigRate, prng, 4.0, maxSpeciesID);
        vector<int> onWhichIsls = archi.findIsl(65);
        assert(onWhichIsls.size() == 2);
        vector<Island> archiCopy = archi.returnArchi();
        int pos = archiCopy[onWhichIsls[0]].findPos(65);
        assert(pos == 1);
        pos = archiCopy[onWhichIsls[1]].findPos(65);
        assert(pos == 0);
        archi.printArchi();
    }
    {
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        SpeciesID maxSpeciesID(n_mainland);
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        const double iniMigRate = pars[1];
        vector<int> happening { 0, 29, 0 };
        archi.doNextEvent(happening, iniMigRate, prng, 4.0, maxSpeciesID);
        vector<int> happening2 { 0, 65, 1 };
        archi.doNextEvent(happening2, iniMigRate, prng, 3.9, maxSpeciesID);
        vector<int> happening3 { 1, 65, 1 };
        archi.doNextEvent(happening3, iniMigRate, prng, 3.8, maxSpeciesID);
        vector<double> sumGloLoc = archi.calculateAllRates(pars, n_mainland, n_islands);
        vector<int> happening4 = archi.sampleNextEvent(sumGloLoc, prng, n_mainland);
        archi.doNextEvent(happening4, iniMigRate, prng, 3.7, maxSpeciesID);
    }
    {
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        SpeciesID maxSpeciesID(n_mainland);
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        const double iniMigRate = pars[1];
        vector<int> happening { 0, 29, 0 };
        archi.doNextEvent(happening, iniMigRate, prng, 4.0, maxSpeciesID);
        vector<int> happening2 { 0, 65, 1 };
        archi.doNextEvent(happening2, iniMigRate, prng, 3.9, maxSpeciesID);
        vector<int> happening3 { 1, 65, 1 };
        archi.doNextEvent(happening3, iniMigRate, prng, 3.8, maxSpeciesID);
        vector<double> sumGloLoc { 1.0, 0.0 };
        vector<int> happening4 = archi.sampleNextEvent(sumGloLoc, prng, n_mainland);
        assert(happening4.size() == 2);  // -> global
        archi.doNextEvent(happening4, iniMigRate, prng, 3.7, maxSpeciesID);
    }
    {  // global clado
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        SpeciesID maxSpeciesID(n_mainland);
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        const double iniMigRate = pars[1];
        vector<int> happening2 { 0, 65, 1 };
        archi.doNextEvent(happening2, iniMigRate, prng, 3.9, maxSpeciesID);
        vector<int> happening3 { 1, 65, 1 };
        archi.doNextEvent(happening3, iniMigRate, prng, 3.8, maxSpeciesID);
        vector<int> happening4 { 5, 65 };
        archi.doNextEvent(happening4, iniMigRate, prng, 3.5, maxSpeciesID);
    }
    {  // global ana
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        SpeciesID maxSpeciesID(n_mainland);
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        const double iniMigRate = pars[1];
        vector<int> happening2 { 0, 65, 1 };
        archi.doNextEvent(happening2, iniMigRate, prng, 3.9, maxSpeciesID);
        vector<int> happening3 { 1, 65, 1 };
        archi.doNextEvent(happening3, iniMigRate, prng, 3.8, maxSpeciesID);
        vector<int> happening4 { 6, 65 };
        archi.doNextEvent(happening4, iniMigRate, prng, 3.5, maxSpeciesID);
    }
    {  // global extinct
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        SpeciesID maxSpeciesID(n_mainland);
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        const double iniMigRate = pars[1];
        vector<int> happening2 { 0, 65, 1 };
        archi.doNextEvent(happening2, iniMigRate, prng, 3.9, maxSpeciesID);
        vector<int> happening3 { 1, 65, 1 };
        archi.doNextEvent(happening3, iniMigRate, prng, 3.8, maxSpeciesID);
        vector<int> happening4 { 7, 65 };
        archi.doNextEvent(happening4, iniMigRate, prng, 3.5, maxSpeciesID);
    }
    {  // local events
        const int n_islands = 2;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        SpeciesID maxSpeciesID(n_mainland);
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        const double iniMigRate = pars[1];
        vector<int> happening2 { 0, 65, 1 };
        archi.doNextEvent(happening2, iniMigRate, prng, 3.9, maxSpeciesID);
        vector<int> happening3 { 1, 65, 1 };
        archi.doNextEvent(happening3, iniMigRate, prng, 3.8, maxSpeciesID);
        vector<int> happening4 { 2, 65, 0 };
        archi.doNextEvent(happening4, iniMigRate, prng, 3.5, maxSpeciesID);
        vector<int> happening5 { 3, 65, 1 };
        archi.doNextEvent(happening5, iniMigRate, prng, 3.4, maxSpeciesID);
        vector<int> happening6 { 4, maxSpeciesID.getMaxSpeciesID(), 1 };
        archi.doNextEvent(happening6, iniMigRate, prng, 3.3, maxSpeciesID);
    }
    {  // more than 2 islands
        const int n_islands = 5;
        const int archi_carryingCap = 50;
        const int n_mainland = 100;
        SpeciesID maxSpeciesID(n_mainland);
        mt19937_64 prng;
        Archipelago archi(n_islands, archi_carryingCap);
        vector<double> pars{ 0.1, 0.1, 0.2, 0.12, 0.3, 0.2, 0.1, 0.12 };
        const double iniMigRate = pars[1];
        vector<int> happening2 { 0, 65, 0 };
        archi.doNextEvent(happening2, iniMigRate, prng, 3.9, maxSpeciesID);
        vector<int> happening3 { 1, 65, 0 };
        archi.doNextEvent(happening3, iniMigRate, prng, 3.8, maxSpeciesID);
        vector<int> happening4 { 0, 65, 2 };
        archi.doNextEvent(happening4, iniMigRate, prng, 3.5, maxSpeciesID);
        vector<int> happening5 { 0, 65, 3 };
        archi.doNextEvent(happening5, iniMigRate, prng, 3.4, maxSpeciesID);
        vector<int> happening6 { 0, 65, 4 };
        archi.doNextEvent(happening6, iniMigRate, prng, 3.3, maxSpeciesID);
    }
}
