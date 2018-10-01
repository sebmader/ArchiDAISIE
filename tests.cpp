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
    { // printing species works
        Species sp1 = Species(0.0, SpeciesID(), SpeciesID(),'M');
        sp1.printSpec();
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
        int destination = island1.drawMigDestinationIsland(0,logGrowthTerms, 0.1, mt19937_64());
        assert(destination == 1);
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
        assert(getEventInt(event) >= 0);
    }
}
