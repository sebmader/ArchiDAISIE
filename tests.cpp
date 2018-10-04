//
// Created by Sebastian Mader on 26.09.2018.
//

#include "tests.h"

using namespace std;

void test_speciesID()
{
    { // default constructor creates speciesID of 0
        SpeciesID spID = SpeciesID();
        assert(spID.getSpeciesID()==0);
    }
    { // incrementing functions is adding 1 to the speciesID
        SpeciesID spID = SpeciesID();
        assert(spID.getSpeciesID()==0);
        spID.incrementSpeciesID();
        assert(spID.getSpeciesID()==1);
    }
    { // creating new speciesID adds 1 to the speciesID and returns copy
        SpeciesID spID = SpeciesID();
        SpeciesID ID = spID.createNewSpeciesID();
        assert(spID.getSpeciesID() == 1);
        assert(ID.getSpeciesID() == 1);
    }
    { // equality operator ==
        SpeciesID spId1 = SpeciesID();
        SpeciesID spId2 = SpeciesID();
        assert(spId1==spId2);
    }
    { // equality operators !=
        SpeciesID spId1 = SpeciesID(1);
        SpeciesID spId2 = SpeciesID(2);
        assert(spId1 != spId2);
    }
    { // equality operators <
        SpeciesID spId1 = SpeciesID(1);
        SpeciesID spId2 = SpeciesID(2);
        assert(spId1<spId2);
    }
    { // equality operators >
        SpeciesID spId1 = SpeciesID(1);
        SpeciesID spId2 = SpeciesID(2);
        assert(spId2>spId1);
    }
    { // equality operators <=
        SpeciesID spId1 = SpeciesID(1);
        SpeciesID spId2 = SpeciesID(2);
        assert(spId1<=spId2);
        spId1.incrementSpeciesID();
        assert(spId1<=spId2);
    }
    { // equality operators >=
        SpeciesID spId1 = SpeciesID(1);
        SpeciesID spId2 = SpeciesID(2);
        assert(spId2>=spId1);
        spId1.incrementSpeciesID();
        assert(spId2>=spId1);
    }
}

void test_species()
{
    { // defauft constructor creates empty species
        Species sp1 = Species();
        assert(sp1.getBirth() == 0.0);
        assert(sp1.getSpecID() == SpeciesID());
        assert(sp1.getParID() == SpeciesID());
        assert(sp1.getStatus() == '0');
    }
    { // setBirth function does set the birth time
        Species sp1 = Species();
        sp1.setBirth(4.0);
        assert(sp1.getBirth() == 4.0);
    }
    { // setStatus function does set the status of species
        Species sp1 = Species();
        sp1.setStatus('I');
        assert(sp1.getStatus() == 'I');
    }
    { // identifies immigrant correctly
        Species sp1 = Species(0.0, SpeciesID(), SpeciesID(),'I');
        assert(sp1.isImmigrant());
    }
    { // identifies migrant correctly
        Species sp1 = Species(0.0, SpeciesID(), SpeciesID(),'M');
        assert(sp1.isMigrant());
    }
}

void test_island()
{
    {   // A default constructed island has a carrying capacity of zero
        const Island island;
        assert(island.getCarryingCap() == 0);
    }
    {   // Setting and getting carrying capacity matches
        const int k{ 12 };
        const Island island(k);
        assert(k==island.getCarryingCap());
    }
    {   // Species cannot be found on empty island
        Island island(10);
        assert(!island.hasSpecies(SpeciesID(42)));
    }
    {   // Adding a species increased the number of species
        Island island(10);
        assert(island.getNSpecies()==0);
        island.addSpecies(Species());
        assert(island.getNSpecies()==1);
    }
    {   // Added species will be found by findSpeciesID
        Island island(10);
        island.addSpecies(Species());
        assert(island.getSpeciesIDs()[0] == SpeciesID());
    }
    {   // The immigration of an absent species increases the number of species
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.getNSpecies()==1);
    }
    {   // The immigration of an absent species adds it to the species vector
        Island island(10);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.hasSpecies(SpeciesID(42)));
    }
    {   // The immigration of an absent species is recognized as an immigrant
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.findSpecies(SpeciesID(42)).getStatus() == 'I');
    }
    {   // An immigrant gets its immigration time stored
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.findSpecies(SpeciesID(42)).getBirth()==3.14);
    }
    {   // When a same species immigrates twice, the last immigration time is stored
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 6.28);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.getNSpecies()==1);
        assert(island.findSpecies(SpeciesID(42)).getBirth()==3.14);
    }
    {  // Species cannot immigrate if island which has as much species as its carrying capacity
        Island island(0);
        try
        {
            island.immigrate(SpeciesID(42), 3.14);
            assert(!"Should not get here"); //!OCLINT accepted idiom
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "Immigration would make number of species"
                                            " exceed carrying capacity.\n");
        }
    }
    {   // Species cannot be found on island before immigration
        Island island(10);
        assert(!island.hasSpecies(SpeciesID(42)));
        island.immigrate(SpeciesID(42), 6.28);
        assert(island.hasSpecies(SpeciesID(42)));
    }
    {   // Extinction decreases the number of species
        Island island(10);
        island.immigrate(SpeciesID(42), 6.28);
        assert(island.getNSpecies()==1);
        island.goExtinct(SpeciesID(42));
        assert(island.getNSpecies()==0);
    }
    { // Extinction of absent species throws an exception
        Island island(1);
        try
        {
            island.goExtinct(SpeciesID(42));
            assert(!"Should not get here"); //!OCLINT accepted idiom
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "Species does not exist on island.\n");
        }
    }
    {   // cladogenesis adds species to island
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        assert(island.getNSpecies() == 1);
        island.speciateClado(SpeciesID(42), 4.0,maxSpeciesID);
        assert(island.getNSpecies() == 2);
    }
    {   // cladogenesis increases maxSpeciesID by 2
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0,maxSpeciesID);
        assert(maxSpeciesID.getSpeciesID() == 50 + 2);
    }
    {   // cladogenesis deletes ancestor
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0,maxSpeciesID);
        assert(!island.hasSpecies(SpeciesID(42)));
    }
    {   // cladogenetic species are of status 'C'
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0,maxSpeciesID);
        assert(island.getSpecies()[0].getStatus() == 'C');
        assert(island.getSpecies()[1].getStatus() == 'C');
    }
    {  // cladogenesis throws exception if species doesn't exists
        Island island(1);
        SpeciesID maxSpeciesID(50);
        try
        {
            island.speciateClado(SpeciesID(42), 4.0,maxSpeciesID);;
            assert(!"Should not get here"); //!OCLINT accepted idiom
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "Species does not exist on island.\n");
        }
    }
    {  // cladogenesis throws exception if it exceeds the carrying cap
        Island island(1);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 4.0);
        try
        {
            island.speciateClado(SpeciesID(42), 4.0,maxSpeciesID);;
            assert(!"Should not get here"); //!OCLINT accepted idiom
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "Cladogenesis would make number of species"
                                            " exceed carrying capacity.\n");
        }
    }
    {   // anagenesis does not add species to island
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        assert(island.getNSpecies() == 1);
        island.speciateAna(SpeciesID(42),maxSpeciesID);
        assert(island.getNSpecies() == 1);
    }
    {   // anagenesis increases maxSpeciesID by 1
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateAna(SpeciesID(42),maxSpeciesID);
        assert(maxSpeciesID.getSpeciesID() == 50 + 1);
    }
    {   // anagenesis deletes ancestor
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateAna(SpeciesID(42),maxSpeciesID);
        assert(!island.hasSpecies(SpeciesID(42)));
    }
    {   // anagenetic species keeps ancestors birth time
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateAna(SpeciesID(42),maxSpeciesID);
        assert(island.findSpecies(maxSpeciesID).getBirth() == 6.28);
    }
    {   // anagenetic species are of status 'A'
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateAna(SpeciesID(42),maxSpeciesID);
        assert(island.getSpecies()[0].getStatus() == 'A');
    }
    {  // anagenesis throws exception if species doesn't exists
        Island island(1);
        SpeciesID maxSpeciesID(50);
        try
        {
            island.speciateAna(SpeciesID(42), maxSpeciesID);;
            assert(!"Should not get here"); //!OCLINT accepted idiom
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "Species does not exist on island.\n");
        }
    }
    {   // migration increases number of species
        Island island(1);
        assert(island.getNSpecies() == 0);
        island.migrate(Species(), 4.0);
        assert(island.getNSpecies() == 1);
    }
    {   // migration makes migration time birth time of migrant
        Island island(1);
        island.migrate(Species(), 4.0);
        assert(island.findSpecies(SpeciesID()).getBirth() == 4.0);
    }
    {   // re-migration of already present species overwrites birth time
        Island island(1);
        island.migrate(Species(), 4.0);
        assert(island.findSpecies(SpeciesID()).getBirth() == 4.0);
        island.migrate(Species(), 2.0);
        assert(island.findSpecies(SpeciesID()).getBirth() == 2.0);
    }
    {   // migration changes status to 'M'
        Island island(1);
        island.migrate(Species(), 4.0);
        assert(island.findSpecies(SpeciesID()).getStatus() == 'M');
    }
    {  // Species cannot migrate if island which has as much species as its carrying capacity
        Island island(0);
        try
        {
            island.migrate(Species(), 3.14);
            assert(!"Should not get here"); //!OCLINT accepted idiom
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "Migration would make number of species"
                                            " exceed carrying capacity.\n");
        }
    }
    {   // drawn island of migration mus be other island than origin
        Island island1(10);
        Island island2(10);
        island1.addSpecies(Species());
        vector<double> logGrowthTerms { getLogGrowth(island1), getLogGrowth(island2)};
        mt19937_64 prng;
        int destination = island1.drawMigDestinationIsland(0,logGrowthTerms, 0.1, prng);
        assert(destination == 1);
    }
    {   // caculating rates of empty island initialises localRates vector
        Island island1(10);
        int n_mainlandSp = 10;
        assert(island1.getLocalRates().size() != 5);
        vector<double> islPars = { 0.1, 0.5, 0.1, 0.05, 0.1 };
        island1.calculateIslRates(islPars, n_mainlandSp, 1, 0);
        assert(island1.getLocalRates().size() == 5);
    }
    {   // caculating rates of empty island leads to immigration
        Island island1(10);
        int n_mainlandSp = 10;
        mt19937_64 prng;
        vector<double> islPars = { 0.1, 0.5, 0.1, 0.05, 0.1 };
        island1.calculateIslRates(islPars, n_mainlandSp, 1, 0);
        event_type event = island1.sampleLocalEvent(prng);
        assert(is_local(event));
        assert(getEventInt(event) == 0);
    }
    {   // add island to another: one immigrant on both islands
        // -> most recent birth time
        Island island1(10);
        Island island2(10);
        island1.immigrate(SpeciesID(), 4.0);
        island2.immigrate(SpeciesID(), 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 2.0);
    }
    {   // add island to another: one immigrant on both islands (other way around)
        // -> most recent birth time
        Island island1(10);
        Island island2(10);
        island1.immigrate(SpeciesID(), 2.0);
        island2.immigrate(SpeciesID(), 4.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 2.0);
    }
    {   // add island to another: immigrant migrated to other island
        // -> colonisation == birth time
        Island island1(10);
        Island island2(10);
        island1.immigrate(SpeciesID(), 4.0);
        island2.migrate(island1.findSpecies(SpeciesID()), 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 4.0);
    }
    {   // add island to another: immigrant migrated to other island (other way around)
        // -> colonisation == birth time
        Island island1(10);
        Island island2(10);
        island2.immigrate(SpeciesID(), 4.0);
        island1.migrate(island2.findSpecies(SpeciesID()), 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 4.0);
    }
    {   // add island to another: cladogenetic sp. migrated to other island
        // -> cladogenesis == birth time
        Island island1(10);
        Island island2(10);
        int n_mainlandSp = 10;
        SpeciesID maxSpeciesID(n_mainlandSp);
        island1.immigrate(SpeciesID(), 4.0);
        island1.speciateClado(SpeciesID(), 3.9, maxSpeciesID);
        island2.migrate(island1.findSpecies(maxSpeciesID), 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 2);
        assert(island1.findSpecies(maxSpeciesID).getBirth() == 3.9);
    }
    {   // add island to another: cladogenetic sp. migrated to other island (other way around)
        // -> cladogenesis == birth time
        Island island1(10);
        Island island2(10);
        int n_mainlandSp = 10;
        SpeciesID maxSpeciesID(n_mainlandSp);
        island2.immigrate(SpeciesID(), 4.0);
        island2.speciateClado(SpeciesID(), 3.9, maxSpeciesID);
        island1.migrate(island2.findSpecies(maxSpeciesID), 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 2);
        assert(island1.findSpecies(maxSpeciesID).getBirth() == 3.9);
    }
    {   // add island to another: anagenetic sp. migrated to other island
        // -> birth time of parent == birth time
        Island island1(10);
        Island island2(10);
        int n_mainlandSp = 10;
        SpeciesID maxSpeciesID(n_mainlandSp);
        island1.immigrate(SpeciesID(), 4.0);
        island1.speciateAna(SpeciesID(), maxSpeciesID);
        island2.migrate(island1.findSpecies(maxSpeciesID), 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(maxSpeciesID).getBirth() == 4.0);
    }
    {   // add island to another: cladogenetic sp. migrated to other island (other way around)
        // -> cladogenesis == birth time
        Island island1(10);
        Island island2(10);
        int n_mainlandSp = 10;
        SpeciesID maxSpeciesID(n_mainlandSp);
        island2.immigrate(SpeciesID(), 4.0);
        island2.speciateAna(SpeciesID(), maxSpeciesID);
        island1.migrate(island2.findSpecies(maxSpeciesID), 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(maxSpeciesID).getBirth() == 4.0);
    }
    {
        Island island(1);
        island.addSpecies(Species());
        island.printIsland();
    }
}


void test_archi()
{
    {  // default constructor creates empty archipelago with 0 islands,
        // 0 carrying cap and no species
        Archipelago archi = Archipelago();
        assert(archi.getNIslands() == 0);
        assert(archi.getNSpecies() == 0);
        assert(archi.getCarryingCap() == 0);
    }
    {  // non-default constructor with carryingCap of 0 and 1 island creates that archipelago
        int n_islands = 1;
        int islCarryingCap = 0;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNIslands() == 1);
        assert(archi.getNSpecies() == 0);
        assert(archi.getCarryingCap() == n_islands*islCarryingCap);
    }
    {  // non-default constructor with 0 islands and a carryingcap creates that archipelago
        int n_islands = 0;
        int islCarryingCap = 1;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNIslands() == 0);
        assert(archi.getNSpecies() == 0);
        assert(archi.getCarryingCap() == n_islands*islCarryingCap);
    }
    {  // non-default constructor creates wanted archipelago
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNIslands() == 2);
        assert(archi.getNSpecies() == 0);
        assert(archi.getCarryingCap() == n_islands*islCarryingCap);
    }
    {  // immigration increases number of species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpecies() == 0);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
    }
    {  // immigration of same species to second islands doesn't increase number of species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpecies() == 0);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                1,
                0.3);
        assert(archi.getNSpecies() == 1);
    }
    {  // after immigration and migration of same species to 2 islands globalSpecies == 1
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getGlobalSpeciesIDs().empty());
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getGlobalSpeciesIDs().empty());
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getGlobalSpeciesIDs().size() == 1);
    }
    {  // after immigration and migration of same species to 2 islands it is found on both islands
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.findIsl(SpeciesID(1)).empty());
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.findIsl(SpeciesID(1)).size() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        vector<int> onWhichIslands = archi.findIsl(SpeciesID(1));
        assert(onWhichIslands.size() == 2);
        assert(onWhichIslands[0] == 0);
        assert(onWhichIslands[1] == 1);
    }
    {  // local cladogenesis increases archi species by one
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpecies() == 0);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 2);
    }
    {  // local anagenesis does not increase archi species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpecies() == 0);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_anagenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
    }
    {  // local extinction decreases archi species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpecies() == 0);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 0);
    }
    {  // doLocalEvent throws if event is not local
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            assert(archi.getNSpecies() == 0);
            int n_mainlandSp = 5;
            SpeciesID maxSpeciesID(n_mainlandSp);
            mt19937_64 prng;
            archi.doLocalEvent(event_type::global_cladogenesis,
                    SpeciesID(1),
                    prng,
                    4.0,
                    maxSpeciesID,
                    0,
                    0.3);
        }
        catch (const exception &e)
        {
            assert(string(e.what()) == "Event is not local.\n");
        }
    }
    // local events with doNextEvent function
    {  // local cladogenesis increases archi species by one
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpecies() == 0);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        archi.doNextEvent(event_type::local_immigration,
                0.3,
                prng,
                4.0,
                maxSpeciesID,
                n_mainlandSp);
        assert(archi.getNSpecies() == 1);
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        archi.doNextEvent(event_type::local_cladogenesis,
                0.3,
                prng,
                3.9,
                maxSpeciesID,
                n_mainlandSp);
        assert(archi.getNSpecies() == 2);
    }
    {  // local anagenesis does not increase archi species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpecies() == 0);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        archi.doNextEvent(event_type::local_immigration,
                0.3,
                prng,
                4.0,
                maxSpeciesID,
                n_mainlandSp);
        assert(archi.getNSpecies() == 1);
        archi.calculateAllRates(iniPars,n_mainlandSp, n_islands);
        archi.doNextEvent(event_type::local_anagenesis,
                0.3,
                prng,
                3.9,
                maxSpeciesID,
                n_mainlandSp);
        assert(archi.getNSpecies() == 1);
    }
    {  // local extinction decreases archi species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpecies() == 0);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        archi.doNextEvent(event_type::local_immigration,
                0.3,
                prng,
                4.0,
                maxSpeciesID,
                n_mainlandSp);
        assert(archi.getNSpecies() == 1);
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        archi.doNextEvent(event_type::local_extinction,
                0.3,
                prng,
                3.9,
                maxSpeciesID,
                n_mainlandSp);
        assert(archi.getNSpecies() == 0);
    }

    {  // calculating rates initialises rates vectors
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        int n_mainlandSp = 5;
        assert(archi.getGlobalRates().empty());
        vector<Island> tmpIsls1 = archi.getIslands();
        assert(tmpIsls1[0].getLocalRates().empty());
        assert(tmpIsls1[1].getLocalRates().empty());
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        assert(archi.getGlobalRates().size() == 3);
        vector<Island> tmpIsls2 = archi.getIslands();
        assert(tmpIsls2[0].getLocalRates().size() == 5);
        assert(tmpIsls2[1].getLocalRates().size() == 5);
    }
    {  // calculated rates are all 0 except for immigration
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        int n_mainlandSp = 5;
        assert(archi.getGlobalRates().empty());
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        assert(archi.getGlobalRates().size() == 3);
        assert(archi.getGlobalRates()[0] == 0.0);
        assert(archi.getGlobalRates()[1] == 0.0);
        assert(archi.getGlobalRates()[2] == 0.0);
        vector<Island> tmpIsls = archi.getIslands();
        assert(tmpIsls[0].getLocalRates()[0] > 0);
        assert(tmpIsls[1].getLocalRates()[0] > 0);
        double sumRatesIsl1 = 0.0;
        double sumRatesIsl2 = 0.0;
        for (int i = 1; i < static_cast<int>(tmpIsls[0].getLocalRates().size()); ++i) {
            sumRatesIsl1 += tmpIsls[0].getLocalRates()[i];
            sumRatesIsl2 += tmpIsls[1].getLocalRates()[i];
        }
        assert(sumRatesIsl1 == 0.0);
        assert(sumRatesIsl2 == 0.0);
    }
    {  // after calculating rates sampled event == immigration
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        int n_mainlandSp = 5;
        assert(archi.getGlobalRates().empty());
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        mt19937_64 prng;
        event_type event = archi.sampleNextEvent(prng);
        assert(event == event_type::local_immigration);
    }
    {   // when there is at least one global species that is immigrant (anagenesis),
        // the global rates are higher than 0
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getGlobalSpeciesIDs().empty());
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getGlobalSpeciesIDs().size() == 1);
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        for (int i = 0; i < (int)archi.getGlobalRates().size(); ++i) {
            assert(archi.getGlobalRates()[i] > 0);
        }
    }
    {   // global cladogenesis increases number of archipelago species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        assert(archi.getNSpecies() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                3.8,
                maxSpeciesID);
        assert(archi.getNSpecies() == 2);
    }
    {   // global cladogenesis doesn't increase the number of species on each island
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        assert(archi.getIslands()[0].getNSpecies() == 0);
        assert(archi.getIslands()[1].getNSpecies() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].getNSpecies() == 1);
        assert(archi.getIslands()[1].getNSpecies() == 0);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].getNSpecies() == 1);
        assert(archi.getIslands()[1].getNSpecies() == 1);
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                3.8,
                maxSpeciesID);
        assert(archi.getIslands()[0].getNSpecies() == 1);
        assert(archi.getIslands()[1].getNSpecies() == 1);
    }
    {   // global cladogenesis creates 2 different species with status 'C'
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].getSpecies()[0].getStatus() == 'I');
        assert(archi.getIslands()[1].getSpecies()[0].getStatus() == 'M');
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                3.8,
                maxSpeciesID);
        assert(archi.getIslands()[0].getSpecies()[0].getStatus() == 'C');
        assert(archi.getIslands()[1].getSpecies()[0].getStatus() == 'C');
        assert(archi.getIslands()[0].getSpecies()[0].getSpecID() !=
                archi.getIslands()[1].getSpecies()[0].getSpecID());
    }
    {   // global cladogenesis throws if species is not present on at least 2 islands
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            int n_mainlandSp = 5;
            SpeciesID maxSpeciesID(n_mainlandSp);
            mt19937_64 prng;
            archi.doGlobalEvent(event_type::global_cladogenesis,
                    SpeciesID(1),
                    prng,
                    3.8,
                    maxSpeciesID);
            assert(!"should not get here!\n"); //!OCLINT
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "Drawn species is present on less than 2 islands. "
                                            "Something's wrong.. (global cladogenesis)\n");
        }
    }
    {   // global anagenesis doesn't increase number of archipelago species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        assert(archi.getNSpecies() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doGlobalEvent(event_type::global_anagenesis,
                SpeciesID(1),
                prng,
                3.8,
                maxSpeciesID);
        assert(archi.getNSpecies() == 1);
    }
    {   // global anagenesis creates one new species with status 'A'
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].getSpecies()[0].getStatus() == 'I');
        assert(archi.getIslands()[1].getSpecies()[0].getStatus() == 'M');
        archi.doGlobalEvent(event_type::global_anagenesis,
                SpeciesID(1),
                prng,
                3.8,
                maxSpeciesID);
        assert(archi.getIslands()[0].getSpecies()[0].getStatus() == 'A');
        assert(archi.getIslands()[1].getSpecies()[0].getStatus() == 'A');
        assert(archi.getIslands()[0].getSpecies()[0].getSpecID()
                == archi.getIslands()[1].getSpecies()[0].getSpecID());
    }
    {   // global anagenesis throws if species is not present on at least 2 islands
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            int n_mainlandSp = 5;
            SpeciesID maxSpeciesID(n_mainlandSp);
            mt19937_64 prng;
            archi.doGlobalEvent(event_type::global_anagenesis,
                    SpeciesID(1),
                    prng,
                    3.8,
                    maxSpeciesID);
            assert(!"should not get here!\n"); //!OCLINT
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "Drawn species is present on less than 2 islands. "
                                            "Something's wrong.. (global anagenesis)\n");
        }
    }
    {   // global extinction decrease number of archipelago species by one
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        assert(archi.getNSpecies() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doGlobalEvent(event_type::global_extinction,
                SpeciesID(1),
                prng,
                3.8,
                maxSpeciesID);
        assert(archi.getNSpecies() == 0);
    }
    {   // global extinction throws exception if species is not on at least 2 islands
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            int n_mainlandSp = 5;
            SpeciesID maxSpeciesID(n_mainlandSp);
            mt19937_64 prng;
            archi.doGlobalEvent(event_type::global_extinction,
                    SpeciesID(1),
                    prng,
                    3.8,
                    maxSpeciesID);
        }
        catch (const exception &e)
        {
            assert(string(e.what()) == "Drawn species is present on less than 2 islands. "
                                        "Something's wrong.. (global extinction)\n");
        }
    }
    {   // doGlobalEvent throws exception if not global event
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            int n_mainlandSp = 5;
            mt19937_64 prng;
            SpeciesID maxSpeciesID(n_mainlandSp);
            archi.doGlobalEvent(event_type::local_immigration,
                    SpeciesID(1),
                    prng,
                    3.8,
                    maxSpeciesID);
            assert(!"Should not get here.\n");  //!OCLINT
        }
        catch (const exception &e)
        {
            assert(string(e.what()) == "Event is not global.\n");
        }
    }
// same with doNextEvent instead of the more specific doLocal / doGlobalEvent
    {   // global cladogenesis increases number of archipelago species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        assert(archi.getNSpecies() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doNextEvent(event_type::global_cladogenesis,
                0.3,
                prng,
                3.8,
                maxSpeciesID, n_mainlandSp);
        assert(archi.getNSpecies() == 2);
    }
    {   // global cladogenesis doesn't increase the number of species on each island
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        assert(archi.getIslands()[0].getNSpecies() == 0);
        assert(archi.getIslands()[1].getNSpecies() == 0);
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].getNSpecies() == 1);
        assert(archi.getIslands()[1].getNSpecies() == 0);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].getNSpecies() == 1);
        assert(archi.getIslands()[1].getNSpecies() == 1);
        archi.doNextEvent(event_type::global_cladogenesis,
                0.3,
                prng,
                3.8,
                maxSpeciesID, n_mainlandSp);
        assert(archi.getIslands()[0].getNSpecies() == 1);
        assert(archi.getIslands()[1].getNSpecies() == 1);
    }
    {   // global cladogenesis creates 2 different species with status 'C'
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].getSpecies()[0].getStatus() == 'I');
        assert(archi.getIslands()[1].getSpecies()[0].getStatus() == 'M');
        archi.doNextEvent(event_type::global_cladogenesis,
                0.3,
                prng,
                3.8,
                maxSpeciesID, n_mainlandSp);
        assert(archi.getIslands()[0].getSpecies()[0].getStatus() == 'C');
        assert(archi.getIslands()[1].getSpecies()[0].getStatus() == 'C');
        assert(archi.getIslands()[0].getSpecies()[0].getSpecID() !=
                archi.getIslands()[1].getSpecies()[0].getSpecID());
    }
    {   // doNextEvent throws if there are no global species (global clado)
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            int n_mainlandSp = 5;
            mt19937_64 prng;
            SpeciesID maxSpeciesID(n_mainlandSp);
            archi.doNextEvent(event_type::global_cladogenesis,
                    0.3,
                    prng,
                    3.8,
                    maxSpeciesID, n_mainlandSp);
            assert(!"should not get here!\n"); //!OCLINT
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "No global species exist on archipelago "
                                            "but drawn event is global.\n");
        }
    }
    {   // global anagenesis doesn't increase number of archipelago species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        assert(archi.getNSpecies() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doNextEvent(event_type::global_anagenesis,
                0.3,
                prng,
                3.8,
                maxSpeciesID, n_mainlandSp);
        assert(archi.getNSpecies() == 1);
    }
    {   // global anagenesis creates one new species with status 'A'
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        mt19937_64 prng;
        SpeciesID maxSpeciesID(n_mainlandSp);
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].getSpecies()[0].getStatus() == 'I');
        assert(archi.getIslands()[1].getSpecies()[0].getStatus() == 'M');
        archi.doNextEvent(event_type::global_anagenesis,
                0.3,
                prng,
                3.8,
                maxSpeciesID, n_mainlandSp);
        assert(archi.getIslands()[0].getSpecies()[0].getStatus() == 'A');
        assert(archi.getIslands()[1].getSpecies()[0].getStatus() == 'A');
        assert(archi.getIslands()[0].getSpecies()[0].getSpecID()
                == archi.getIslands()[1].getSpecies()[0].getSpecID());
    }
    {   // doNextEvent throws if there are no global species (global ana)
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            int n_mainlandSp = 5;
            SpeciesID maxSpeciesID(n_mainlandSp);
            mt19937_64 prng;
            archi.doNextEvent(event_type::global_anagenesis,
                    0.3,
                    prng,
                    3.8,
                    maxSpeciesID, n_mainlandSp);
            assert(!"should not get here!\n"); //!OCLINT
        }
        catch (const std::exception& e)
        {
            assert(std::string(e.what()) == "No global species exist on archipelago "
                                            "but drawn event is global.\n");
        }
    }
    {   // global extinction decrease number of archipelago species by one
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        assert(archi.getNSpecies() == 0);
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        archi.doNextEvent(event_type::global_extinction,
                0.3,
                prng,
                3.8,
                maxSpeciesID, n_mainlandSp);
        assert(archi.getNSpecies() == 0);
    }
    {   // doNextEvent throws if there are no global species (global extinct)
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            int n_mainlandSp = 5;
            mt19937_64 prng;
            SpeciesID maxSpeciesID(n_mainlandSp);
            archi.doNextEvent(event_type::global_extinction,
                    0.3,
                    prng,
                    3.8,
                    maxSpeciesID, n_mainlandSp);
        }
        catch (const exception &e)
        {
            assert(string(e.what()) == "No global species exist on archipelago "
                                       "but drawn event is global.\n");
        }
    }
    {   // first event = immigration with doNextEvent function
        // increases number of species on archipelago
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        mt19937_64 prng;
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        event_type event = archi.sampleNextEvent(prng);
        assert(archi.getNSpecies() == 0);
        archi.doNextEvent(event, iniPars[1], prng,
                4.0, maxSpeciesID, n_mainlandSp);
        assert(archi.getNSpecies() == 1);
    }
    {   // archipelago with 3 islands, species inhabiting all,
        // doesn't create duplicates in global species vector
        int n_islands = 3;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        mt19937_64 prng;
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        event_type event = archi.sampleNextEvent(prng);
        assert(archi.getNSpecies() == 0);
        archi.doNextEvent(event, iniPars[1], prng,
                4.0, maxSpeciesID, n_mainlandSp);
        assert(archi.getNSpecies() == 1);
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        assert(archi.getGlobalSpeciesIDs().empty());
        archi.doNextEvent(event_type::local_migration,
                0.3,
                prng,
                3.5,
                maxSpeciesID,
                n_mainlandSp);
        assert(archi.getGlobalSpeciesIDs().size() == 1);
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        archi.doNextEvent(event_type::local_migration,
                0.3,
                prng,
                3.3,
                maxSpeciesID,
                n_mainlandSp);
        assert(archi.getNSpecies() == 1);
        assert(archi.getGlobalSpeciesIDs().size() == 1);
        archi.calculateAllRates(iniPars, n_mainlandSp, n_islands);
        archi.doNextEvent(event_type::local_migration,
                0.3,
                prng,
                3.2,
                maxSpeciesID,
                n_mainlandSp);
    }
}
