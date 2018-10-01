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
        int destination = island1.drawMigDestinationIsland(0,logGrowthTerms, 0.1, mt19937_64());
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
    {  // default constructor creates empty archipelago with 0 islands
        // and 0 carrying cap
        Archipelago archi = Archipelago();
        assert(archi.getNIslands() == 0);
        assert(archi.getCarryingCap() == 0);
    }
    {  //

    }
}
