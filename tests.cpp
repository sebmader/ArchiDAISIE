//
// Created by Sebastian Mader on 26.09.2018.
//

#include "tests.h"

using namespace std;


void test_other_functions() //!OCLINT indeed long function, don't care it is a test
{
    // extractSumOfRates:
    { // extracting the sum of rates of empty island
        Island island(10);
        assert(extractSumOfRates(island)==0);
    }
    { // extracting the sum of rates of island with calculated rates
        Island island(10);
        const int n_mainlandSp = 5;
        vector<double> pars({ 0.1,0.1,0.1,0.2,0.1 });
        island.calculateIslRates(pars,n_mainlandSp,1,0.1);
        assert(extractSumOfRates(island)==n_mainlandSp*pars[0]);
    }
    // getLogGrowth
    { // empty island
        Island island(10);
        assert(getLogGrowth(island)==1);
    }
    { // colonised island
        Island island(10);
        island.immigrate(SpeciesID(),4.0);
        assert(getLogGrowth(island)==0.9);
    }
    // whichMainAncestor
    { // doesn't find any on emtpy island
        Island island = Island();
        assert(whichMainAncestors(island.getSpecies()).empty());
    }
    { // does find 1 on island with 1 species
        Island island = Island(10);
        island.addSpecies(Species());
        assert(whichMainAncestors(island.getSpecies()).size()==1);
    }
    { // just finds 1 on island with 2 species that are sisters
        Island island = Island(10);
        Species sp1 = Species(3.2,SpeciesID(3),SpeciesID(102),'0',false,3.2,3.2);
        Species sp2 = Species(3.0,SpeciesID(3),SpeciesID(103),'0',false,3.2,3.2);
        island.addSpecies(sp1);
        island.addSpecies(sp2);
        assert(whichMainAncestors(island.getSpecies()).size()==1);
        assert(whichMainAncestors(island.getSpecies())[0].isSister(sp1));
        assert(whichMainAncestors(island.getSpecies())[0].isSister(sp2));
        assert(whichMainAncestors(island.getSpecies())[0] != sp1);
    }
    { // if 2 species from same mainland ancestor but with different colonisation times
        // then the older colonisation time will be saved
        Island island = Island(10);
        Species sp1 = Species(1.0,SpeciesID(3),SpeciesID(102),'0',false,3.2,3.2);
        Species sp2 = Species(1.0,SpeciesID(3),SpeciesID(103),'0',false,1.0,1.0);
        island.addSpecies(sp1);
        island.addSpecies(sp2);
        assert(whichMainAncestors(island.getSpecies()).size()==1);
        assert(whichMainAncestors(island.getSpecies())[0].isSister(sp1));
        assert(!whichMainAncestors(island.getSpecies())[0].isSister(sp2));
        assert(whichMainAncestors(island.getSpecies())[0] != sp1);
    }
    { // if 3 species from same mainland ancestor but with different colonisation times
        // then the older colonisation time will be saved
        Island island = Island(10);
        Species sp1 = Species(1.0,SpeciesID(3),SpeciesID(102),'0',false,3.2,3.2);
        Species sp2 = Species(1.0,SpeciesID(3),SpeciesID(103),'0',false,1.0,1.0);
        Species sp3 = Species(1.0,SpeciesID(3),SpeciesID(106),'0',false,2.3,2.3);
        island.addSpecies(sp2);
        island.addSpecies(sp3);
        island.addSpecies(sp1);
        assert(whichMainAncestors(island.getSpecies()).size()==1);
        assert(whichMainAncestors(island.getSpecies())[0].isSister(sp1));
        assert(!whichMainAncestors(island.getSpecies())[0].isSister(sp2));
        assert(!whichMainAncestors(island.getSpecies())[0].isSister(sp3));
        assert(whichMainAncestors(island.getSpecies())[0] != sp1);
    }
    // outputBranching
    /*
    { // output of branching considers re-immigrant of existing clade as missing species
        // and counts it as missing species -> works!
        Island island = Island(10);
        const int n_mainland = 20;
        SpeciesID maxSpeciesID = SpeciesID(n_mainland);
        island.immigrate(SpeciesID(),2.0);
        island.speciateClado(SpeciesID(),1.5,maxSpeciesID);
        island.immigrate(SpeciesID(),1.0);
        ofstream ofs("test.txt");
        outputBranching(island,ofs);
    }
    */
}

void test_speciesID() //!OCLINT indeed long function, don't care it is a test
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

void test_species() //!OCLINT indeed long function, don't care it is a test
{
    { // defauft constructor creates empty species
        Species sp1 = Species();
        assert(sp1.getBirth() == 0.0);
        assert(sp1.getSpecID() == SpeciesID());
        assert(sp1.getParID() == SpeciesID());
        assert(sp1.getStatus() == '0');
        assert(!sp1.hasMigrated());
        assert(sp1.getAncestralBT() == 0.0);
        assert(sp1.getColonisationT() == 0.0);
        assert(sp1.getCladoStates().empty());
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
    { // setAncestralBT function does set the birth time of the clade
        Species sp1 = Species();
        sp1.setAncestralBT(2.0);
        assert(sp1.getAncestralBT() == 2.0);
    }
    { // identifies immigrant correctly
        Species sp1 = Species();
        sp1.setStatus('I');
        assert(sp1.isImmigrant());
    }
    { // identifies cladogenetic species correctly -> only about cladoStates, not status
        Species sp1 = Species(3.5,SpeciesID(),SpeciesID(5),'A',false,3.5,3.5,{'a','b'});
        assert(sp1.isCladogenetic());
    }
    { // operator== identifies equal species
        Species sp1 = Species();
        Species sp2 = Species();
        assert(sp1 == sp2);
        sp1 = Species(3.5,SpeciesID(),SpeciesID(5),'A',false,3.5,3.5,{'a','b'});
        sp2 = Species(3.5,SpeciesID(),SpeciesID(5),'A',false,3.5,3.5,{'a','b'});
        assert(sp1 == sp2);
    }
    { // operator!= identifies unequal species
        Species sp1 = Species();
        Species sp2 = Species(1.0);
        assert(sp1 != sp2);
        sp1 = Species(3.5,SpeciesID(),SpeciesID(5),'A',false,3.5,3.5,{'a','b'});
        sp2 = Species(3.5,SpeciesID(),SpeciesID(6),'A',false,3.5,3.5,{'a','b'});
        assert(sp1 != sp2);
        sp1 = Species(3.5,SpeciesID(),SpeciesID(5),'A',false,3.5,3.5,{'a','b'});
        sp2 = Species(3.5,SpeciesID(),SpeciesID(5),'A',false,3.5,3.5,{'b','a'});
        assert(sp1 != sp2);
    }
    { // species is not recognized as its own sister
        Species sp1 = Species(3.2, SpeciesID(5), SpeciesID(53), 'A', false,
                3.2, 3.2, {'a'});
        assert(!sp1.isSister(sp1));
    }
    { // identifies sister species correctly
        Species sp1 = Species(3.2, SpeciesID(5), SpeciesID(53), 'A', false,
                3.2, 3.2, {'a'});
        Species sp2 = Species(2.6, SpeciesID(5), SpeciesID(52), 'C', true,
                3.4, 3.2, {'b', 'a'});
        assert(sp1.isSister(sp2));
    }
}

void test_island() //!OCLINT indeed long function, don't care it is a test
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
    {   // A constructed island has 0 colonisations
        const Island island;
        assert(island.getNColonisations() == 0);
    }
    {   // Species cannot be found on empty island
        Island island(10);
        assert(!island.hasSpecies(SpeciesID()));
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
    {   // An absent species that immigrates is recognized as an immigrant
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.findSpecies(SpeciesID(42)).isImmigrant());
    }
    {   // An immigrant gets its immigration time stored as birth time
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
    {   // An immigrant gets its immigration time stored as clade birth time
        Island island(10);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.findSpecies(SpeciesID(42)).getAncestralBT()==3.14);
    }
    {   // An immigrant gets its immigration time stored as colonisation time
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.findSpecies(SpeciesID(42)).getColonisationT()==3.14);
    }
    {   // An immigrant has not migrated
        Island island(10);
        assert(island.getNSpecies()==0);
        island.immigrate(SpeciesID(42), 3.14);
        assert(!island.findSpecies(SpeciesID(42)).hasMigrated());
    }
    {   // An immigrant has empty cladogenesis states
        Island island(10);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.findSpecies(SpeciesID(42)).getCladoStates().empty());
    }
    {   // Immigration increases colonisations of island
        Island island(10);
        assert(island.getNColonisations() == 0);
        island.immigrate(SpeciesID(42), 3.14);
        assert(island.getNColonisations() == 1);
    }
    {   // Species cannot be found on island before immigration
        Island island(10);
        assert(!island.hasSpecies(SpeciesID(42)));
        island.immigrate(SpeciesID(42), 6.28);
        assert(island.hasSpecies(SpeciesID(42)));
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
    {   // 1st daughter inherits birth time of parent
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        assert(island.findSpecies(SpeciesID(51)).getBirth() == 6.28);
    }
    {   // 2st daughter gets new birth time
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        assert(island.findSpecies(SpeciesID(52)).getBirth() == 4.0);
    }
    {   // cladogenesis saves the clado states (size == 1)
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        assert(island.findSpecies(SpeciesID(51)).getCladoStates().size() == 1);
        assert(island.findSpecies(SpeciesID(52)).getCladoStates().size() == 1);
    }
    {   // 1st daughter gets cladogenesis state 'a'
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        vector<char> cladoStates = island.findSpecies(SpeciesID(51)).getCladoStates();
        assert(cladoStates[0] == 'a');
    }
    {   // 2st daughter gets cladogenesis state 'b'
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        vector<char> cladoStates = island.findSpecies(SpeciesID(52)).getCladoStates();
        assert(cladoStates[0] == 'b');
    }
    {   // cladogenetic daughters are recognized as sisters
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        assert(island.findSpecies(SpeciesID(51)).isSister(island.findSpecies(SpeciesID(52))));
    }
    {   // clades from different immigrations from same mainland species
        // are not recognized as sisters
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        island.immigrate(SpeciesID(42), 6.0);
        assert(!island.findSpecies(SpeciesID(51)).isSister(island.findSpecies(SpeciesID(42))));
        assert(!island.findSpecies(SpeciesID(52)).isSister(island.findSpecies(SpeciesID(42))));
    }
    {   // double cladogenenetic species are all sisters
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        island.speciateClado(SpeciesID(51), 3.7, maxSpeciesID);
        Species sp1 = island.findSpecies(SpeciesID(52));
        Species sp2 = island.findSpecies(SpeciesID(53));
        Species sp3 = island.findSpecies(SpeciesID(54));
        assert(sp1.isSister(sp2) && sp1.isSister(sp3));
        assert(sp2.isSister(sp1) && sp2.isSister(sp3));
        assert(sp3.isSister(sp1) && sp3.isSister(sp2));
    }
    {   // double cladogenesis: identifies most recent sister species/clade correctly
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        island.speciateClado(SpeciesID(51), 3.7, maxSpeciesID);
        Species sp1 = island.findSpecies(SpeciesID(52));
        Species sp2 = island.findSpecies(SpeciesID(53));
        Species sp3 = island.findSpecies(SpeciesID(54));
        assert(sp1.isMostRecentSis(sp2) && sp1.isMostRecentSis(sp3));
        assert(sp2.isMostRecentSis(sp3) && sp3.isMostRecentSis(sp2));
        assert(!sp2.isMostRecentSis(sp1) && !sp3.isMostRecentSis(sp1));
    }
    {   // double cladogenesis saves both clado states
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        island.speciateClado(SpeciesID(51), 4.0, maxSpeciesID);
        assert(island.findSpecies(SpeciesID(53)).getCladoStates().size() == 2);
        assert(island.findSpecies(SpeciesID(54)).getCladoStates().size() == 2);
    }
    {   // double clado: 1st daughter of 1st daughter gets cladogenesis states 'a','a'
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        island.speciateClado(SpeciesID(51), 4.0, maxSpeciesID);
        vector<char> cladoStates = island.findSpecies(SpeciesID(53)).getCladoStates();
        assert(cladoStates[0] == 'a' && cladoStates[1] == 'a');
    }
    {   // double clado: 2st daughter of 1st daughter gets cladogenesis states 'a','b'
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        island.speciateClado(SpeciesID(51), 4.0, maxSpeciesID);
        vector<char> cladoStates = island.findSpecies(SpeciesID(54)).getCladoStates();
        assert(cladoStates[0] == 'a' && cladoStates[1] == 'b');
    }
    {   // double clado: 1st daughter of 2nd daughter gets cladogenesis states 'b','a'
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        island.speciateClado(SpeciesID(52), 4.0, maxSpeciesID);
        vector<char> cladoStates = island.findSpecies(SpeciesID(53)).getCladoStates();
        assert(cladoStates[0] == 'b' && cladoStates[1] == 'a');
    }
    {   // double clado: 2st daughter of 2nd daughter gets cladogenesis states 'b','b'
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        island.speciateClado(SpeciesID(52), 4.0, maxSpeciesID);
        vector<char> cladoStates = island.findSpecies(SpeciesID(54)).getCladoStates();
        assert(cladoStates[0] == 'b' && cladoStates[1] == 'b');
    }
    {   // 1st daughter inherits ancestral birth
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        assert(island.findSpecies(SpeciesID(51)).getAncestralBT() == 6.28);
    }
    {   // 2nd daughter gets new ancestral birth = its own birth time
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0, maxSpeciesID);
        assert(island.findSpecies(SpeciesID(52)).getAncestralBT() == 4.0);
        assert(island.findSpecies(SpeciesID(52)).getBirth()
                == island.findSpecies(SpeciesID(52)).getAncestralBT());
    }
    {   // cladogenetic species are of status 'C'
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0,maxSpeciesID);
        assert(island.getSpecies()[0].getStatus() == 'C');
        assert(island.getSpecies()[1].getStatus() == 'C');
    }
    {   // cladogenetic species both have the old species' colonisation time
        Island island(10);
        SpeciesID maxSpeciesID(50);
        island.immigrate(SpeciesID(42), 6.28);
        island.speciateClado(SpeciesID(42), 4.0,maxSpeciesID);
        assert(island.getSpecies()[0].getColonisationT() == 6.28);
        assert(island.getSpecies()[1].getColonisationT() == 6.28);
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
    {   // anagenetic species saves ancestors birth time as ancestral BT
        Island island(10);
        SpeciesID maxSpeciesID(50);
        Species sp1 = Species(1.0,SpeciesID(),SpeciesID(5),'I',false, 2.0,3.0,{'a'});
        island.addSpecies(sp1);
        island.speciateAna(SpeciesID(5),maxSpeciesID);
        assert(island.findSpecies(maxSpeciesID).getAncestralBT() == sp1.getBirth());
    }
    {   // anagenetic species keeps ancestors colonisation time
        Island island(10);
        SpeciesID maxSpeciesID(50);
        Species sp1 = Species(1.0,SpeciesID(),SpeciesID(5),'I',false, 2.0,3.0,{'a'});
        island.addSpecies(sp1);
        island.speciateAna(SpeciesID(5),maxSpeciesID);
        assert(island.findSpecies(maxSpeciesID).getColonisationT() == sp1.getColonisationT());
    }
    {   // anagenetic species keeps ancestors cladogenesis states
        Island island(10);
        SpeciesID maxSpeciesID(50);
        Species sp1 = Species(1.0, SpeciesID(), SpeciesID(), '0', false, 1.0, 1.0,
                { 'a','b','a' });
        island.addSpecies(sp1);
        island.speciateAna(SpeciesID(),maxSpeciesID);
        assert(island.findSpecies(maxSpeciesID).getCladoStates().size() == 3);
        assert(island.findSpecies(maxSpeciesID).getCladoStates() == sp1.getCladoStates());
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
    {   // migration of absent species saves old birthT as ancestral birthT
        Island island1(1);
        Island island2(1);
        island1.immigrate(SpeciesID(), 4.0);
        island2.migrate(island1.findSpecies(SpeciesID()), 2.0);
        assert(island2.findSpecies(SpeciesID()).getAncestralBT()
                == island1.findSpecies(SpeciesID()).getBirth());
    }
    {   // migration of absent species keeps colonisation time
        Island island1(1);
        Island island2(1);
        island1.immigrate(SpeciesID(), 4.0);
        island2.migrate(island1.findSpecies(SpeciesID()), 2.0);
        assert(island2.findSpecies(SpeciesID()).getColonisationT()==4.0);
    }
    {   // migration of absent species adds new daughter state 'b' to migrant
        Island island1(1);
        island1.immigrate(SpeciesID(),3.0);
        Island island2(1);
        island2.migrate(island1.findSpecies(SpeciesID()), 2.0);
        assert(island2.findSpecies(SpeciesID()).getCladoStates().size() == 1);
        assert(island2.findSpecies(SpeciesID()).getCladoStates()[0] == 'b');
    }
    {   // migration does not change status
        Island island(1);
        island.migrate(Species(), 4.0);
        assert(island.findSpecies(SpeciesID()).getStatus() == '0');
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
    {   // drawn island of migration must be other island than origin
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
    {   // add island to another: one pop of mainland species on both islands
        // -> most recent coloniser stays (-> re-immigration)
        Island island1(10);
        Island island2(10);
        island1.immigrate(SpeciesID(), 4.0);
        island2.immigrate(SpeciesID(), 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 2.0);
    }
    {   // add island to another: one pop of mainland species on both islands (other way around)
        // -> most recent coloniser stays (-> re-immigration)
        Island island1(10);
        Island island2(10);
        island1.immigrate(SpeciesID(), 2.0);
        island2.immigrate(SpeciesID(), 4.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 2.0);
    }
    {   // add 3 islands together: one pop of mainland species on both islands
        // -> oldest pop of most recent coloniser stays (-> re-immigration)
        Island island1(10);
        Island island2(10);
        Island island3(10);
        island1.immigrate(SpeciesID(), 4.0);
        island2.immigrate(SpeciesID(), 2.0);
        island3.migrate(island2.findSpecies(SpeciesID()),1.0);
        island1.addIsland(island3);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 1.0);
        assert(island1.findSpecies(SpeciesID()).getColonisationT() == 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 2.0);
        assert(island1.findSpecies(SpeciesID()).getColonisationT() == 2.0);
    }
    {   // add 3 islands together: one pop of mainland species on both islands (other way around)
        // -> oldest pop of most recent coloniser stays (-> re-immigration)
        Island island1(10);
        Island island2(10);
        Island island3(10);
        island3.immigrate(SpeciesID(), 4.0);
        island2.immigrate(SpeciesID(), 2.0);
        island1.migrate(island2.findSpecies(SpeciesID()),1.0);
        island1.addIsland(island3);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 1.0);
        assert(island1.findSpecies(SpeciesID()).getColonisationT() == 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 2.0);
        assert(island1.findSpecies(SpeciesID()).getColonisationT() == 2.0);
    }
    {   // add 3 islands together: independent pops of mainland species on each islands
        // -> oldest pop of most recent coloniser stays (-> re-immigration)
        Island island1(10);
        Island island2(10);
        Island island3(10);
        island1.immigrate(SpeciesID(), 4.0);
        island2.immigrate(SpeciesID(), 2.0);
        island3.migrate(island2.findSpecies(SpeciesID()),1.0);
        island2.immigrate(SpeciesID(), 0.9);
        island1.addIsland(island3);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 1.0);
        assert(island1.findSpecies(SpeciesID()).getColonisationT() == 2.0);
        island1.addIsland(island2);
        assert(island1.getNSpecies() == 1);
        assert(island1.findSpecies(SpeciesID()).getBirth() == 0.9);
        assert(island1.findSpecies(SpeciesID()).getColonisationT() == 0.9);
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
    {   // add island to another: anagenetic sp. migrated to other island (other way around)
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
}


void test_archi() //!OCLINT indeed long function, don't care it is a test
{
    {  // default constructor creates empty archipelago with 0 islands,
        // 0 carrying cap, no species and 0 colonisations
        Archipelago archi = Archipelago();
        assert(archi.getNIslands() == 0);
        assert(archi.getNSpeciesID() == 0);
        assert(archi.getCarryingCap() == 0);
        assert(archi.getNColonisations() == 0);
    }
    {  // non-default constructor with carryingCap of 0 and 1 island creates that archipelago

        int n_islands = 1;
        int islCarryingCap = 0;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNIslands() == 1);
        assert(archi.getNSpeciesID() == 0);
        assert(archi.getCarryingCap() == n_islands*islCarryingCap);
        assert(archi.getNColonisations() == 0);
    }
    {  // non-default constructor with 0 islands and a carryingcap creates that archipelago
        int n_islands = 0;
        int islCarryingCap = 1;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNIslands() == 0);
        assert(archi.getNSpeciesID() == 0);
        assert(archi.getCarryingCap() == n_islands*islCarryingCap);
    }
    {  // non-default constructor creates wanted archipelago
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNIslands() == 2);
        assert(archi.getNSpeciesID() == 0);
        assert(archi.getCarryingCap() == n_islands*islCarryingCap);
    }
    {  // immigration increases number of species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpeciesID() == 0);
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
    {  // immigration increases number of colonisations
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNColonisations() == 0);
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
        assert(archi.getNColonisations() == 1);
    }
    {  // immigration of same species to second islands doesn't increase number of speciesIDs
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpeciesID() == 0);
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
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                1,
                0.3);
        assert(archi.getNSpeciesID() == 1);
    }
    {  // re-immigration of same species replaces it
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpeciesID() == 0);
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
                3.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getBirth() == 3.0);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getColonisationT() == 3.0);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getAncestralBT() == 3.0);
    }
    {  // re-immigration of same species replaces migrant and corrects cladoStates
        // of population on original island of origin
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpeciesID() == 0);
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
                3.5,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates()[0]=='b');
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                3.0,
                maxSpeciesID,
                1,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().empty());
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().empty());
    }
    {  // re-immigration of same species replaces migrant and corrects cladoStates
        // of population on original island of origin (more cladoStates)
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpeciesID() == 0);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        Species sp1 = Species(4.0,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{ 'b','b','b' });
        archi.addSpecies(sp1,0);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.5,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size()==4);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[3]=='a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().size()==4);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates()[3]=='b');
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                3.0,
                maxSpeciesID,
                1,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size()==3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[2]=='b');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().empty());
    }
    {  // migration does increase number of species
        // but not SpeciesIDs
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
        assert(archi.getNSpecies() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.5,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpecies() == 2);
        assert(archi.getNSpeciesID() == 1);
    }
    {   // re-migration of already present species overwrites birth time
        // if migrant is older than resident
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
                3.5,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[1].hasSpecies(SpeciesID(1)));
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getBirth() == 3.0);
    }
    {   // re-migration of already present species doesn't overwrite birth time
        // if migrant is younger than resident
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
                3.5,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.0,
                maxSpeciesID,
                1,
                0.3);
        assert(archi.getIslands()[0].hasSpecies(SpeciesID(1)));
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getBirth() == 4.0);
    }
    {  // migration does add clado state 'a' to population on island of origin
        // and 'b' to population on island of destination
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
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().empty());
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.5,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size() == 1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[0] == 'a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().size() == 1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates()[0] == 'b');
    }
    {   // re-migration of already present species does not add a new daughter state
        // if migrant comes from the same origin
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
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().empty());
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.5,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size() == 1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[0] == 'a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().size() == 1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates()[0] == 'b');
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size() == 1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[0] == 'a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().size() == 1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates()[0] == 'b');
    }
    {   // re-migration of already present species does correct daughter state of original
        // migrant if new migrant (re-migrant) comes from another island of origin
        // TODO: how to test that?
    }
    {  // migration does replace already present species if migrant is re-immigrant
        // = more recent colonisation
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
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                3.0,
                maxSpeciesID,
                1,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                2.5,
                maxSpeciesID,
                1,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getBirth() == 2.5);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getColonisationT() == 3.0);
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
        assert(archi.getNSpeciesID() == 0);
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
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 2);
    }
    {  // local anagenesis does not increase archi species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpeciesID() == 0);
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
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_anagenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
    }
    {  // local extinction decreases archi species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpeciesID() == 0);
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
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 0);
    }
    {  // doLocalEvent throws if event is not local
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            assert(archi.getNSpeciesID() == 0);
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
        assert(archi.getNSpeciesID() == 0);
        int n_mainlandSp(5);
        vector<SpeciesID> mainlandSpIDs(n_mainlandSp);
        for(size_t i = 0; i < mainlandSpIDs.size(); ++i) {
            mainlandSpIDs[i] = SpeciesID(i);
        }
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        archi.doNextEvent(event_type::local_immigration,
                0.3,
                prng,
                4.0,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        archi.doNextEvent(event_type::local_cladogenesis,
                0.3,
                prng,
                3.9,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 2);
    }
    {  // local anagenesis does not increase archi species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpeciesID() == 0);
        int n_mainlandSp = 5;
        vector<SpeciesID> mainlandSpIDs(n_mainlandSp);
        for(size_t i = 0; i < mainlandSpIDs.size(); ++i) {
            mainlandSpIDs[i] = SpeciesID(i);
        }
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        archi.doNextEvent(event_type::local_immigration,
                0.3,
                prng,
                4.0,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
        archi.calculateAllRates(iniPars,mainlandSpIDs, n_islands);
        archi.doNextEvent(event_type::local_anagenesis,
                0.3,
                prng,
                3.9,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
    }
    {  // local extinction decreases archi species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        assert(archi.getNSpeciesID() == 0);
        int n_mainlandSp = 5;
        vector<SpeciesID> mainlandSpIDs(n_mainlandSp);
        for(size_t i = 0; i < mainlandSpIDs.size(); ++i) {
            mainlandSpIDs[i] = SpeciesID(i);
        }
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        archi.doNextEvent(event_type::local_immigration,
                0.3,
                prng,
                4.0,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        archi.doNextEvent(event_type::local_extinction,
                0.3,
                prng,
                3.9,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 0);
    }
    {  // calculating rates initialises rates vectors
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        int n_mainlandSp = 5;
        vector<SpeciesID> mainlandSpIDs(n_mainlandSp);
        for(size_t i = 0; i < mainlandSpIDs.size(); ++i) {
            mainlandSpIDs[i] = SpeciesID(i);
        }
        assert(archi.getGlobalRates().empty());
        vector<Island> tmpIsls1 = archi.getIslands();
        assert(tmpIsls1[0].getLocalRates().empty());
        assert(tmpIsls1[1].getLocalRates().empty());
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
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
        vector<SpeciesID> mainlandSpIDs(n_mainlandSp);
        for(size_t i = 0; i < mainlandSpIDs.size(); ++i) {
            mainlandSpIDs[i] = SpeciesID(i);
        }
        assert(archi.getGlobalRates().empty());
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
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
        vector<SpeciesID> mainlandSpIDs(n_mainlandSp);
        for(size_t i = 0; i < mainlandSpIDs.size(); ++i) {
            mainlandSpIDs[i] = SpeciesID(i);
        }
        assert(archi.getGlobalRates().empty());
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
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
        vector<SpeciesID> mainlandSpIDs(n_mainlandSp);
        for(size_t i = 0; i < mainlandSpIDs.size(); ++i) {
            mainlandSpIDs[i] = SpeciesID(i);
        }
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
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
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
        assert(archi.getNSpeciesID() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        assert(archi.getNSpeciesID() == 2);
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
                maxSpeciesID);
        assert(archi.getIslands()[0].getNSpecies() == 1);
        assert(archi.getIslands()[1].getNSpecies() == 1);
    }
    // TODO: test for birth times and ancestral birth times
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
        assert(archi.getIslands()[0].getSpecies()[0].isImmigrant());
        assert(archi.getIslands()[1].getSpecies()[0].hasMigrated());
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        assert(archi.getIslands()[0].getSpecies()[0].isCladogenetic());
        assert(archi.getIslands()[1].getSpecies()[0].isCladogenetic());
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
            archi.addSpecies(Species(1.0,SpeciesID(),SpeciesID(1)),1);
            archi.speciateGlobalClado(SpeciesID(1),prng,maxSpeciesID);
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
        assert(archi.getNSpeciesID() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doGlobalEvent(event_type::global_anagenesis,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        assert(archi.getNSpeciesID() == 1);
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
        assert(archi.getIslands()[0].getSpecies()[0].isImmigrant());
        assert(archi.getIslands()[1].getSpecies()[0].hasMigrated());
        archi.doGlobalEvent(event_type::global_anagenesis,
                SpeciesID(1),
                prng,
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
            archi.addSpecies(Species(1.0, SpeciesID(),SpeciesID(1)),0);
            archi.speciateGlobalAna(SpeciesID(1),maxSpeciesID);
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
        assert(archi.getNSpeciesID() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        assert(archi.isGlobal(SpeciesID(1)));
        archi.doGlobalEvent(event_type::global_extinction,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        assert(archi.getNSpeciesID() == 0);
    }
    {   // global extinction throws exception if species is not on at least 2 islands
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            int n_mainlandSp = 5;
            SpeciesID maxSpeciesID(n_mainlandSp);
            mt19937_64 prng;
            archi.addSpecies(Species(2.0,SpeciesID(),SpeciesID(1)),0);
            archi.goGlobalExtinct(SpeciesID(1));
        }
        catch (const exception &e)
        {
            assert(string(e.what()) == "Drawn species is present on less than 2 islands. "
                                       "Something's wrong.. (global extinction)\n");
        }
    }
    {   // doGlobalEvent throws exception if species is not on at least 2 islands
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
                    maxSpeciesID);
        }
        catch (const exception &e)
        {
            assert(string(e.what()) == "Drawn species is present on less than 2 islands. "
                                        "Something's wrong.. (doGlobalEvent)\n");
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
            archi.addSpecies(Species(1.1,SpeciesID(),SpeciesID(1)),0);
            archi.addSpecies(Species(1.0,SpeciesID(),SpeciesID(1)),1);
            archi.doGlobalEvent(event_type::local_immigration,
                    SpeciesID(1),
                    prng,
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
        assert(archi.getNSpeciesID() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doNextEvent(event_type::global_cladogenesis,
                0.3,
                prng,
                3.8,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 2);
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
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
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
        assert(archi.getIslands()[0].getSpecies()[0].isImmigrant());
        assert(archi.getIslands()[1].getSpecies()[0].hasMigrated());
        archi.doNextEvent(event_type::global_cladogenesis,
                0.3,
                prng,
                3.8,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getIslands()[0].getSpecies()[0].isCladogenetic());
        assert(archi.getIslands()[1].getSpecies()[0].isCladogenetic());
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
                    maxSpeciesID,
                    vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
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
        assert(archi.getNSpeciesID() == 0);
        mt19937_64 prng;
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doNextEvent(event_type::global_anagenesis,
                0.3,
                prng,
                3.8,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
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
        assert(archi.getIslands()[0].getSpecies()[0].isImmigrant());
        assert(archi.getIslands()[1].getSpecies()[0].hasMigrated());
        archi.doNextEvent(event_type::global_anagenesis,
                0.3,
                prng,
                3.8,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
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
                    maxSpeciesID,
                    vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
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
        assert(archi.getNSpeciesID() == 0);
        archi.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getNSpeciesID() == 1);
        archi.doNextEvent(event_type::global_extinction,
                0.3,
                prng,
                3.8,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 0);
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
                    maxSpeciesID,
                    vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        }
        catch (const exception &e)
        {
            assert(string(e.what()) == "No global species exist on archipelago "
                                       "but drawn event is global.\n");
        }
    }
    {   // first event = immigration with doNextEvent function
        // increases number of species on archipelago and colonisations
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        vector<SpeciesID> mainlandSpIDs(n_mainlandSp);
        for(size_t i = 0; i < mainlandSpIDs.size(); ++i) {
            mainlandSpIDs[i] = SpeciesID(i);
        }
        SpeciesID maxSpeciesID(n_mainlandSp);
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        mt19937_64 prng;
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        event_type event = archi.sampleNextEvent(prng);
        assert(archi.getNSpeciesID() == 0);
        archi.doNextEvent(event, iniPars[1], prng,
                4.0, maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
        assert(archi.getNColonisations() == 1);
    }
    {   // archipelago with 3 islands, species inhabiting all,
        // doesn't create duplicates in global species vector
        int n_islands = 3;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        vector<SpeciesID> mainlandSpIDs(n_mainlandSp);
        for(size_t i = 0; i < mainlandSpIDs.size(); ++i) {
            mainlandSpIDs[i] = SpeciesID(i);
        }
        SpeciesID maxSpeciesID(n_mainlandSp);
        vector<double> iniPars { 0.05, 0.5, 0.2, 0.1, 0.2, 0.1, 0.05, 0.1 };
        mt19937_64 prng;
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        event_type event = archi.sampleNextEvent(prng);
        assert(archi.getNSpeciesID() == 0);
        archi.doNextEvent(event, iniPars[1], prng,
                4.0, maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
        assert(archi.getGlobalSpeciesIDs().empty());
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        archi.doNextEvent(event_type::local_migration,
                0.3,
                prng,
                3.5,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
        assert(archi.getGlobalSpeciesIDs().size() == 1);
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        archi.doNextEvent(event_type::local_migration,
                0.3,
                prng,
                3.3,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
        assert(archi.getGlobalSpeciesIDs().size() == 1);
        archi.calculateAllRates(iniPars, mainlandSpIDs, n_islands);
        archi.doNextEvent(event_type::local_migration,
                0.3,
                prng,
                3.2,
                maxSpeciesID,
                vector<SpeciesID>( { SpeciesID(n_mainlandSp) } ));
        assert(archi.getNSpeciesID() == 1);
        assert(archi.getGlobalSpeciesIDs().size() == 1);
    }
    // most recent sister function
    {  // after migration migrant is found as sister of ancestral population (and vice versa)
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
                3.9,
                maxSpeciesID,
                0,
                0.3);
        Species pop1 = archi.getIslands()[0].findSpecies(SpeciesID(1));
        Species pop2 = archi.getIslands()[1].findSpecies(SpeciesID(1));
        assert(archi.findMostRecentSistersPops(pop1).size() == 1);
        assert(archi.findMostRecentSistersPops(pop1)[0] == pop2);
        assert(archi.findMostRecentSistersPops(pop2).size() == 1);
        assert(archi.findMostRecentSistersPops(pop2)[0] == pop1);
    }
    {  // after double migration (-> 3 islands) the younger ones still
        // just find each other as most recent sisters, not the older one
        // older finds both younger ones though
        int n_islands = 3;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{'a'}),0);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','a'}),1);
        archi.addSpecies(Species(3.8,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','b'}),2);
        Species pop1 = archi.getIslands()[0].getSpecies()[0];
        Species pop2 = archi.getIslands()[1].getSpecies()[0];
        Species pop3 = archi.getIslands()[2].getSpecies()[0];
        assert(archi.findMostRecentSistersPops(pop1).size() == 2);
        assert(archi.findMostRecentSistersPops(pop1)[0] == pop2);
        assert(archi.findMostRecentSistersPops(pop1)[1] == pop3);
        assert(archi.findMostRecentSistersPops(pop2).size() == 1);
        assert(archi.findMostRecentSistersPops(pop2)[0] == pop3);
        assert(archi.findMostRecentSistersPops(pop3).size() == 1);
        assert(archi.findMostRecentSistersPops(pop3)[0] != pop1);
    }
    {  // after local cladogenesis daughters find each other as most recent sisters
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        Species sp = archi.getIslands()[0].findSpecies(SpeciesID(6));
        Species sp1 = archi.getIslands()[0].findSpecies(SpeciesID(7));
        assert(archi.findMostRecentSistersPops(sp).size() == 1);
        assert(archi.findMostRecentSistersPops(sp)[0] == sp1);
        assert(archi.findMostRecentSistersPops(sp1).size() == 1);
        assert(archi.findMostRecentSistersPops(sp1)[0] == sp);
    }
    {  // after local cladogenesis, migration of daughter 'b'
        // most recent sister populations are found on both islands separately
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(7),
                prng,
                3.8,
                maxSpeciesID,
                0,
                0.3);
        Species sp = archi.getIslands()[0].findSpecies(SpeciesID(6));
        assert(archi.findMostRecentSistersPops(sp).size() == 2);
        assert(archi.findMostRecentSistersPops(sp)[0].getSpecID() == SpeciesID(7));
        assert(archi.findMostRecentSistersPops(sp)[1].getSpecID() == SpeciesID(7));
    }
    {  // after local cladogenesis, migration of daughter 'b'
        // daughter 'b' populations just find each other as sisters, not daughter 'a'
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(7),
                prng,
                3.8,
                maxSpeciesID,
                0,
                0.3);
        const Species spB1 = archi.getIslands()[0].findSpecies(SpeciesID(7));
        const Species spB2 = archi.getIslands()[1].findSpecies(SpeciesID(7));
        assert(archi.findMostRecentSistersPops(spB1).size() == 1);
        assert(archi.findMostRecentSistersPops(spB1)[0] == spB2);
        assert(archi.findMostRecentSistersPops(spB2).size() == 1);
        assert(archi.findMostRecentSistersPops(spB2)[0] == spB1);
    }
    {  // after local cladogenesis, migration of daughter 'b'
        // daughter 'a' population finds both daughter 'b' populations as MRS
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(7),
                prng,
                3.8,
                maxSpeciesID,
                0,
                0.3);
        const Species spB1 = archi.getIslands()[0].findSpecies(SpeciesID(7));
        const Species spB2 = archi.getIslands()[1].findSpecies(SpeciesID(7));
        const Species spA = archi.getIslands()[0].findSpecies(SpeciesID(6));
        assert(archi.findMostRecentSistersPops(spA).size() == 2);
        assert(archi.findMostRecentSistersPops(spA)[0] == spB1);
        assert(archi.findMostRecentSistersPops(spA)[1] == spB2);
    }
    {  // after local cladogenesis, migration & extinction of first daughter 'b'
        // second daughter 'b' is found as most recent sister of daughter 'a' and vice versa
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(7),
                prng,
                3.8,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(7),
                prng,
                3.7,
                maxSpeciesID,
                0,
                0.3);
        Species spA = archi.getIslands()[0].findSpecies(SpeciesID(6));
        Species spB = archi.getIslands()[1].findSpecies(SpeciesID(7));
        assert(archi.findMostRecentSistersPops(spA).size() == 1);
        assert(archi.findMostRecentSistersPops(spA)[0] == spB);
        assert(archi.findMostRecentSistersPops(spB).size() == 1);
        assert(archi.findMostRecentSistersPops(spB)[0] == spA);
    }
    {  // after global cladogenesis daughters find other as most recent sisters
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
                3.8,
                maxSpeciesID,
                0,
                0.3);
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        Species spA = archi.getIslands()[0].findSpecies(SpeciesID(6));
        Species spB = archi.getIslands()[1].findSpecies(SpeciesID(7));
        assert(archi.findMostRecentSistersPops(spA).size() == 1);
        assert(archi.findMostRecentSistersPops(spA)[0] == spB);
        assert(archi.findMostRecentSistersPops(spB).size() == 1);
        assert(archi.findMostRecentSistersPops(spB)[0] == spA);
    }
    {  // after global cladogenesis of 3 daughters (-> 3 islands) the younger ones still
        // just find each other as most recent sisters, not the older one
        int n_islands = 3;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{'a'}),0);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','a'}),1);
        archi.addSpecies(Species(3.8,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','b'}),2);
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        Species spA = archi.getIslands()[0].getSpecies()[0];
        Species spAB = archi.getIslands()[1].getSpecies()[0];
        Species spB = archi.getIslands()[2].getSpecies()[0];
        assert(archi.findMostRecentSistersPops(spA).size() == 2);
        assert(archi.findMostRecentSistersPops(spA)[0] == spAB);
        assert(archi.findMostRecentSistersPops(spA)[1] == spB);
        assert(archi.findMostRecentSistersPops(spAB).size() == 1);
        assert(archi.findMostRecentSistersPops(spAB)[0] == spB);
        assert(archi.findMostRecentSistersPops(spB).size() == 1);
        assert(archi.findMostRecentSistersPops(spB)[0] == spAB);
    }
    // testing correcting sisters after extinction
    { // local migration
        // if older pop of a species dies, the younger one inherits the birth time (BT)
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
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getBirth() == 4.0);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getBirth() == 3.6);
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(1),
                prng,
                3.2,
                maxSpeciesID,
                0,
                0.3);
        assert(!archi.getIslands()[0].hasSpecies(SpeciesID(1)));
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getBirth() == 4.0);
    }
    { // local migration
        // if older pop of a species dies, younger one looses resp. daughter state
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
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates()[0]=='b');
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(1),
                prng,
                3.2,
                maxSpeciesID,
                0,
                0.3);
        assert(!archi.getIslands()[0].hasSpecies(SpeciesID(1)));
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().empty());
    }
    { // local migration
        // if younger pop of a species dies, BT stays the same
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
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getBirth() == 4.0);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getBirth() == 3.6);
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(1),
                prng,
                3.2,
                maxSpeciesID,
                1,
                0.3);
        assert(!archi.getIslands()[1].hasSpecies(SpeciesID(1)));
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getBirth() == 4.0);
    }
    { // local migration
        // if younger pop of a species dies, older looses daughter state
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
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates()[0]=='b');
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(1),
                prng,
                3.2,
                maxSpeciesID,
                1,
                0.3);
        assert(!archi.getIslands()[1].hasSpecies(SpeciesID(1)));
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().empty());
    }
    { // local cladogenesis
        // if older sister species goes extinct, the younger one inherits the birth time (BT)
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getBirth()==4.0);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getBirth()==3.6);
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(6),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getBirth()==4.0);
    }
    { // local cladogenesis
        // if older sister species goes extinct, the younger one looses daughter state
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getCladoStates()[0]=='b');
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(6),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getCladoStates().empty());
    }
    { // local cladogenesis
        // if younger sister species goes extinct, the older one keeps its BT
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getBirth()==4.0);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getBirth()==3.6);
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(7),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getBirth()==4.0);
    }
    { // local cladogenesis
        // if younger sister species goes extinct, the older one looses daughter state
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getCladoStates()[0]=='b');
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(7),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates().empty());
    }
    { // global cladogenesis
        // if older sister of a species dies, the younger one inherits the birth time (BT)
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
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getBirth() == 4.0);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getBirth() == 3.6);
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getBirth() == 4.0);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getBirth() == 3.6);
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(6),
                prng,
                3.2,
                maxSpeciesID,
                0,
                0.3);
        assert(!archi.getIslands()[0].hasSpecies(SpeciesID(6)));
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getBirth() == 4.0);
    }
    { // global cladogenesis
        // if older sister of a species dies, the younger one looses daughter state
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
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().empty());
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates()[0]=='b');
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getCladoStates().size()==1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getCladoStates()[0]=='b');
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(6),
                prng,
                3.2,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getCladoStates().empty());
    }
    { // global cladogenesis
        // if younger sister of a species dies, the older keeps its birth time (BT)
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
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getBirth() == 4.0);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getBirth() == 3.6);
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getBirth() == 4.0);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getBirth() == 3.6);
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(7),
                prng,
                3.2,
                maxSpeciesID,
                1,
                0.3);
        assert(!archi.getIslands()[1].hasSpecies(SpeciesID(7)));
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getBirth() == 4.0);
    }
    { // global cladogenesis
        // if yougner sister of a species dies, the older one looses daughter state
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
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().empty());
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(1),
                prng,
                3.6,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(1)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates().size()==1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(1)).getCladoStates()[0]=='b');
        archi.doGlobalEvent(event_type::global_cladogenesis,
                SpeciesID(1),
                prng,
                maxSpeciesID);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getCladoStates().size()==1);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getCladoStates()[0]=='b');
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(7),
                prng,
                3.2,
                maxSpeciesID,
                1,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates().empty());
    }
    {  // after local cladogenesis, migration & extinction of first daughter 'b'
        // second daughter 'b' looses daughter state
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
        archi.doLocalEvent(event_type::local_cladogenesis,
                SpeciesID(1),
                prng,
                3.9,
                maxSpeciesID,
                0,
                0.3);
        archi.doLocalEvent(event_type::local_migration,
                SpeciesID(7),
                prng,
                3.8,
                maxSpeciesID,
                0,
                0.3);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates().size()==1);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(6)).getCladoStates()[0]=='a');
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getCladoStates().size()==2);
        assert(archi.getIslands()[0].findSpecies(SpeciesID(7)).getCladoStates()
                ==vector<char>({'b','a'}));
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getCladoStates().size()==2);
        assert(archi.getIslands()[1].findSpecies(SpeciesID(7)).getCladoStates()
                ==vector<char>({'b','b'}));
        archi.doLocalEvent(event_type::local_extinction,
                SpeciesID(7),
                prng,
                3.7,
                maxSpeciesID,
                0,
                0.3);
        Species spA = archi.getIslands()[0].findSpecies(SpeciesID(6));
        Species spB = archi.getIslands()[1].findSpecies(SpeciesID(7));
        assert(spA.getCladoStates().size()==1);
        assert(spA.getCladoStates()[0] == 'a');
        assert(spB.getCladoStates().size()==1);
        assert(spB.getCladoStates()[0] == 'b');
    }

    { // getGlobalSpeciesIDs does not create duplicates if present on 2 islands
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{'a'}),0);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','a'}),1);
        assert(archi.getGlobalSpeciesIDs().size()==1);
        assert(archi.getGlobalSpeciesIDs()[0]==SpeciesID(1));
    }
    { // getGlobalSpeciesIDs does not create duplicates if present on 3 islands
        int n_islands = 3;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{'a'}),0);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','a'}),1);
        archi.addSpecies(Species(3.8,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','b'}),2);
        assert(archi.getGlobalSpeciesIDs().size()==1);
        assert(archi.getGlobalSpeciesIDs()[0]==SpeciesID(1));
    }
    { // getGlobalSpeciesIDs does not create duplicates if present on 4 islands
        int n_islands = 4;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{'a'}),0);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','a'}),1);
        archi.addSpecies(Species(3.8,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','b','a'}),2);
        archi.addSpecies(Species(3.5,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b','b','b'}),3);
        assert(archi.getGlobalSpeciesIDs().size()==1);
        assert(archi.getGlobalSpeciesIDs()[0]==SpeciesID(1));
    }
    { // getGlobalSpeciesIDs recognises multiple global species on differnt islands
        int n_islands = 4;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{'a'}),0);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b'}),1);
        archi.addSpecies(Species(3.8,SpeciesID(2),SpeciesID(2),'I',false,4.0,4.0,{'a'}),2);
        archi.addSpecies(Species(3.5,SpeciesID(2),SpeciesID(2),'I',true,4.0,4.0,{'b'}),3);
        assert(archi.getGlobalSpeciesIDs().size()==2);
        assert(archi.getGlobalSpeciesIDs()[0]==SpeciesID(1));
        assert(archi.getGlobalSpeciesIDs()[1]==SpeciesID(2));
    }
    { // getGlobalSpeciesIDs recognises multiple global species on same islands
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{'a'}),0);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{'b'}),1);
        archi.addSpecies(Species(3.8,SpeciesID(2),SpeciesID(2),'I',false,4.0,4.0,{'a'}),0);
        archi.addSpecies(Species(3.5,SpeciesID(2),SpeciesID(2),'I',true,4.0,4.0,{'b'}),1);
        assert(archi.getGlobalSpeciesIDs().size()==2);
        assert(archi.getGlobalSpeciesIDs()[0]==SpeciesID(1));
        assert(archi.getGlobalSpeciesIDs()[1]==SpeciesID(2));
    }
    { // getGlobalSpeciesIDs recognises multiple global endemic species on same islands
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(6),'C',false,4.0,4.0,{'a','a'}),0);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(6),'C',true,4.0,4.0,{'a','b'}),1);
        archi.addSpecies(Species(3.8,SpeciesID(2),SpeciesID(7),'C',false,3.8,3.8,{'b','a'}),0);
        archi.addSpecies(Species(3.5,SpeciesID(2),SpeciesID(7),'C',true,3.8,3.8,{'b','b'}),1);
        assert(archi.getGlobalSpeciesIDs().size()==2);
        assert(archi.getGlobalSpeciesIDs()[0]==SpeciesID(6));
        assert(archi.getGlobalSpeciesIDs()[1]==SpeciesID(7));
    }
    // adding/consolidating archipelagos
    { // summed archipelago has species that both archipelagos had
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi1 = Archipelago(n_islands, islCarryingCap);
        Species sp1 = Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{});
        archi1.addSpecies(sp1,0);
        Archipelago archi2 = Archipelago(n_islands, islCarryingCap);
        Species sp2 = Species(4.0,SpeciesID(2),SpeciesID(2),'I',false,4.0,4.0,{});
        archi2.addSpecies(sp2,1);
        archi1.addArchi(archi2);
        assert(archi1.getNSpeciesID() == 2);
        assert(archi1.hasSpecies(SpeciesID(1)));
        assert(archi1.findIslSpecies(sp1.getSpecID())[0] == sp1);
        assert(archi1.hasSpecies(SpeciesID(2)));
        assert(archi1.findIslSpecies(sp2.getSpecID())[0] == sp2);
    }
    { // summed archipelago has species that both archipelagos had
        // different kinds of species
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi1 = Archipelago(n_islands, islCarryingCap);
        Species sp1 = Species(4.0,SpeciesID(1),SpeciesID(1),'C',false,4.0,4.0,{ 'a' });
        archi1.addSpecies(sp1,0);
        Archipelago archi2 = Archipelago(n_islands, islCarryingCap);
        Species sp2 = Species(3.0,SpeciesID(2),SpeciesID(2),'C',true,4.0,4.0,{ 'b' });
        archi2.addSpecies(sp2,1);
        archi1.addArchi(archi2);
        assert(archi1.getNSpeciesID() == 2);
        assert(archi1.hasSpecies(SpeciesID(1)));
        assert(archi1.findIslSpecies(sp1.getSpecID())[0] == sp1);
        assert(archi1.hasSpecies(SpeciesID(2)));
        assert(archi1.findIslSpecies(sp2.getSpecID())[0] == sp2);
    }
    { // summed archipelago has just one species present when same species is
        // present on the same island on respective archipelago
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi1 = Archipelago(n_islands, islCarryingCap);
        Species sp1 = Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{});
        archi1.addSpecies(sp1, 0);
        Archipelago archi2 = Archipelago(n_islands, islCarryingCap);
        Species sp2 = Species(3.4,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{});
        archi2.addSpecies(sp2, 0);
        archi1.addArchi(archi2);
        assert(archi1.getNSpecies() == 1);
        assert(archi1.hasSpecies(SpeciesID(1)));
        assert(archi1.findIslSpecies(SpeciesID(1))[0] == sp1);
    }
    { // summed archipelago has just one species present when same species is
        // present on the same island on respective archipelago
        // different kind of species -> older BT but younger coloT stays
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi1 = Archipelago(n_islands, islCarryingCap);
        Species sp1 = Species(2.0,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{});
        archi1.addSpecies(sp1, 0);
        Archipelago archi2 = Archipelago(n_islands, islCarryingCap);
        Species sp2 = Species(2.8,SpeciesID(1),SpeciesID(1),'I',true,3.4,3.4,{});
        archi2.addSpecies(sp2, 0);
        archi1.addArchi(archi2);
        assert(archi1.getNSpecies() == 1);
        assert(archi1.hasSpecies(SpeciesID(1)));
        assert(archi1.findIslSpecies(SpeciesID(1))[0] == sp2);
    }
    { // summed archipelago has both species when same species is present
        // on different islands on respective archipelagos
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi1 = Archipelago(n_islands, islCarryingCap);
        archi1.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{}),0);
        Archipelago archi2 = Archipelago(n_islands, islCarryingCap);
        archi2.addSpecies(Species(3.4,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{}),1);
        archi1.addArchi(archi2);
        assert(archi1.getNSpeciesID() == 1);
        assert(archi1.getNSpecies() == 2);
        assert(archi1.hasSpecies(SpeciesID(1)));
    }
    { // consolidated archipelago has species present on single islands
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        Species sp1 = Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{});
        archi.addSpecies(sp1,0);
        Species sp2 = Species(3.5,SpeciesID(2),SpeciesID(2),'I',false,3.5,3.5,{});
        archi.addSpecies(sp2,1);
        Island consolidatedArchi = archi.makeArchiTo1Island();
        assert(consolidatedArchi.getNSpecies()==2);
        assert(consolidatedArchi.hasSpecies(SpeciesID(1)));
        assert(consolidatedArchi.findSpecies(SpeciesID(1)) == sp1);
        assert(consolidatedArchi.hasSpecies(SpeciesID(2)));
        assert(consolidatedArchi.findSpecies(SpeciesID(2)) == sp2);
    }
    { // consolidated archipelago deletes duplicates -> older pop stays
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        Species sp1 = Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{});
        archi.addSpecies(sp1,0);
        Species sp2 = Species(3.5,SpeciesID(1),SpeciesID(1),'I',true,4.0,4.0,{});
        archi.addSpecies(sp2,1);
        Island consolidatedArchi = archi.makeArchiTo1Island();
        assert(consolidatedArchi.getNSpecies()==1);
        assert(consolidatedArchi.hasSpecies(SpeciesID(1)));
        assert(consolidatedArchi.findSpecies(SpeciesID(1)) == sp1);
    }
    { // consolidated archipelago deletes duplicates -> younger coloniser stays
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        Species sp1 = Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{});
        archi.addSpecies(sp1,0);
        Species sp2 = Species(3.5,SpeciesID(1),SpeciesID(1),'I',false,3.5,3.5,{});
        archi.addSpecies(sp2,1);
        Island consolidatedArchi = archi.makeArchiTo1Island();
        assert(consolidatedArchi.getNSpecies()==1);
        assert(consolidatedArchi.hasSpecies(SpeciesID(1)));
        assert(consolidatedArchi.findSpecies(SpeciesID(1)) == sp2);
    }
    { // summed archipelago has sum of colonisations of both archis
        int n_islands = 2;
        int islCarryingCap = 5;
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        Archipelago archi1 = Archipelago(n_islands, islCarryingCap);
        assert(archi1.getNColonisations() == 0);
        archi1.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi1.getNColonisations() == 1);
        Archipelago archi2 = Archipelago(n_islands, islCarryingCap);
        assert(archi2.getNColonisations() == 0);
        archi2.doLocalEvent(event_type::local_immigration,
                SpeciesID(2),
                prng,
                3.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi2.getNColonisations() == 1);
        archi1.addArchi(archi2);
        assert(archi1.getNColonisations() == 2);
    }
    { // summed archipelago has sum of colonisations of both archis
        // also if same species
        int n_islands = 2;
        int islCarryingCap = 5;
        int n_mainlandSp = 5;
        SpeciesID maxSpeciesID(n_mainlandSp);
        mt19937_64 prng;
        Archipelago archi1 = Archipelago(n_islands, islCarryingCap);
        assert(archi1.getNColonisations() == 0);
        archi1.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                4.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi1.getNColonisations() == 1);
        Archipelago archi2 = Archipelago(n_islands, islCarryingCap);
        assert(archi2.getNColonisations() == 0);
        archi2.doLocalEvent(event_type::local_immigration,
                SpeciesID(1),
                prng,
                3.0,
                maxSpeciesID,
                0,
                0.3);
        assert(archi2.getNColonisations() == 1);
        archi1.addArchi(archi2);
        assert(archi1.getNColonisations() == 2);
    }
}

void test_STT()
{
    { // empty constructor creates emtpy STT
        STT stt = STT();
        assert(stt.getTime() == 0.0);
        assert(stt.getNImmigrants() == 0);
        assert(stt.getNAnagenetic() == 0);
        assert(stt.getNCladogenetic() == 0);
        assert(stt.getNColonisations() == 0);
    }
    { // constructor assigns values correctly
        STT stt(2.0, 2, 1, 3, 1);
        assert(stt.getTime() == 2.0);
        assert(stt.getNImmigrants() == 2);
        assert(stt.getNAnagenetic() == 1);
        assert(stt.getNCladogenetic() == 3);
        assert(stt.getNColonisations() == 1);
    }
}

void test_STTtable() //!OCLINT indeed long function, don't care it is a test
{
    { // empty constructor creates empty STTtable
        STTtable sttTable = STTtable();
        assert(sttTable.size() == 0);
    }
    { // constructor creates STTtable with input as size and those rows being 0s
        STTtable sttTable = STTtable(1);
        assert(sttTable.size() == 1);
        assert(sttTable.getSTTtable()[0].getTime() == 0.0);
        assert(sttTable.getSTTtable()[0].getNImmigrants() == 0);
        assert(sttTable.getSTTtable()[0].getNAnagenetic() == 0);
        assert(sttTable.getSTTtable()[0].getNCladogenetic() == 0);
        assert(sttTable.getSTTtable()[0].getNColonisations() == 0);
    }
    { // constructor creates STTtable with input as size
        // and those rows being second input
        STTtable sttTable = STTtable(1, STT(2.0, 2, 1, 3, 1));
        assert(sttTable.size() == 1);
        assert(sttTable.getSTTtable()[0].getTime() == 2.0);
        assert(sttTable.getSTTtable()[0].getNImmigrants() == 2);
        assert(sttTable.getSTTtable()[0].getNAnagenetic() == 1);
        assert(sttTable.getSTTtable()[0].getNCladogenetic() == 3);
        assert(sttTable.getSTTtable()[0].getNColonisations() == 1);
    }
    { // update function: present species increases num of colonisations
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(5.0,SpeciesID(),SpeciesID(),'I'),0);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateSTTtable(archi,4.0);
        assert(sttTable.size()==2);
        assert(sttTable.getSTTtable()[1].getTime()==4.0);
        assert(sttTable.getSTTtable()[1].getNColonisations()==1);
    }
    { // update function: 2 present species increases num of colonisations to 2
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(5.0,SpeciesID(6),SpeciesID(6),'I',false,6.0,6.0),0);
        archi.addSpecies(Species(5.0,SpeciesID(7),SpeciesID(7),'I',false,5.5,5.5),0);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateSTTtable(archi,5.0);
        assert(sttTable.size()==2);
        assert(sttTable.getSTTtable()[1].getTime()==5.0);
        assert(sttTable.getSTTtable()[1].getNColonisations()==2);
    }
    { // update function: sister species are identified as one colonisation
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(5.0,SpeciesID(),SpeciesID(1),'I',false,5.0,5.0),0);
        archi.addSpecies(Species(4.9,SpeciesID(),SpeciesID(2),'I',false,5.0,5.0),1);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateSTTtable(archi,4.0);
        assert(sttTable.size()==2);
        assert(sttTable.getSTTtable()[1].getTime()==4.0);
        assert(sttTable.getSTTtable()[1].getNColonisations()==1);
        // output:
        // cout << sttTable << '\n';
    }
    { // update function identifies and saves immigrant species on archipelago correctly
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,4.0,{}),0);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateSTTtable(archi,4.0);
        assert(sttTable.size()==2);
        assert(sttTable.getSTTtable()[1].getTime() == 4.0);
        assert(sttTable.getSTTtable()[1].getNImmigrants() == 1);
        assert(sttTable.getSTTtable()[1].getNAnagenetic() == 0);
        assert(sttTable.getSTTtable()[1].getNCladogenetic() == 0);
        assert(sttTable.getSTTtable()[1].getNColonisations()==1);
    }
    { // update function identifies and saves anagenetic species on archipelago correctly
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'A',false,4.0,4.0,{}),0);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateSTTtable(archi,4.0);
        assert(sttTable.size()==2);
        assert(sttTable.getSTTtable()[1].getTime() == 4.0);
        assert(sttTable.getSTTtable()[1].getNImmigrants() == 0);
        assert(sttTable.getSTTtable()[1].getNAnagenetic() == 1);
        assert(sttTable.getSTTtable()[1].getNCladogenetic() == 0);
        assert(sttTable.getSTTtable()[1].getNColonisations()==1);
    }
    { // update function identifies and saves cladogenetic species on archipelago correctly
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(4.0,SpeciesID(1),SpeciesID(1),'C',false,4.0,4.0,{}),0);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateSTTtable(archi,4.0);
        assert(sttTable.size()==2);
        assert(sttTable.getSTTtable()[1].getTime() == 4.0);
        assert(sttTable.getSTTtable()[1].getNImmigrants() == 0);
        assert(sttTable.getSTTtable()[1].getNAnagenetic() == 0);
        assert(sttTable.getSTTtable()[1].getNCladogenetic() == 1);
        assert(sttTable.getSTTtable()[1].getNColonisations()==1);
    }
    { // update function throws if status of species is default ('0')
        try {
            int n_islands = 2;
            int islCarryingCap = 5;
            Archipelago archi = Archipelago(n_islands, islCarryingCap);
            archi.addSpecies(Species(),0);
            STTtable sttTable = STTtable(1);
            assert(sttTable.size()==1);
            sttTable.updateSTTtable(archi,4.0);
        }
        catch (exception& ex) {
            assert(string(ex.what()) == "Status of species is unknown.\n");
        }
    }
    { // operator "[]" works for STT class
        STT firstRow = STT(3.4,1,6,3,3);
        STTtable stt = STTtable(1,firstRow);
        assert(stt[0] == firstRow);
    }
    { // const operator "[]" works for STT class
        const STT firstRow = STT(3.4,1,6,3,3);
        const STTtable stt = STTtable(1,firstRow);
        assert(stt[0] == firstRow);
    }
    { // merging of STTtables has size == sample size input + 1
        const STT secondRow = STT(3.4,1,6,3,3);
        STTtable stt = STTtable(2,STT(4.0));
        stt[1] = secondRow;
        const STT secondRow2 = STT(3.0,1,1,2,3);
        STTtable stt2 = STTtable(2,STT(4.0));
        stt2[1] = secondRow2;
        const int n_slices = 8;
        STTtable merge = mergeSTTtables({ stt,stt2 },n_slices);
        assert(merge[0].getTime() == 4.0);
        assert(merge[1].getTime() == 4.0 - (4.0/n_slices));
        assert(merge.getSTTtable().back().getTime()== 0.0);
        assert(merge.size()==n_slices+1);
    }
    { // merging of STTtables results in time steps of islandAge/n_slices
        // and ends at zero (last elements time)
        const STT secondRow = STT(3.4,1,6,3,3);
        STTtable stt = STTtable(2,STT(4.0));
        stt[1] = secondRow;
        const STT secondRow2 = STT(3.0,1,1,2,3);
        STTtable stt2 = STTtable(2,STT(4.0));
        stt2[1] = secondRow2;
        const int n_slices = 8;
        STTtable merge = mergeSTTtables({ stt,stt2 },n_slices);
        assert(merge[0].getTime() == 4.0);
        assert(merge[1].getTime() == 4.0 - (4.0/n_slices));
        assert(merge.getSTTtable().back().getTime()== 0.0);
    }
    { // merging of STTtables results in sums of single STT inputs at specific sample times
        const STT secondRow = STT(3.6,1,6,3,3);
        STTtable stt = STTtable(2,STT(4.0));
        stt[1] = secondRow;
        const STT secondRow2 = STT(3.1,1,1,2,3);
        STTtable stt2 = STTtable(2,STT(4.0));
        stt2[1] = secondRow2;
        const int n_slices = 8;
        STTtable merge = mergeSTTtables({ stt,stt2 },n_slices);
        assert(merge[1].getTime()==3.5);
        assert(merge[1].getNImmigrants() == secondRow.getNImmigrants());
        assert(merge[1].getNAnagenetic() == secondRow.getNAnagenetic());
        assert(merge[1].getNCladogenetic() == secondRow.getNCladogenetic());
        assert(merge[1].getNColonisations() == secondRow.getNColonisations());
        assert(merge[2].getTime()==3.0);
        assert(merge[2].getNImmigrants() == secondRow.getNImmigrants()
                + secondRow2.getNImmigrants());
        assert(merge[2].getNAnagenetic() == secondRow.getNAnagenetic()
                + secondRow2.getNAnagenetic());
        assert(merge[2].getNCladogenetic() == secondRow.getNCladogenetic()
                + secondRow2.getNCladogenetic());
        assert(merge[2].getNColonisations() == secondRow.getNColonisations()
                + secondRow2.getNColonisations());
    }
    // problem: merging of island STT tables counts each island population as a single species
      // no way of identifying duplicate species on basis of STT tables or on final L table
    // possible solution: STT table updated during run-time is archipelago-wide
      // -> use mergingArchi function in update function

    {  // STT update does not delete 'A' duplicates
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'A',false,4.0,4.0,{}),0);
        archi.addSpecies(Species(3.0,SpeciesID(1),SpeciesID(1),'A',true,4.0,4.0,{}),1);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateSTTtable(archi,3.0);
        assert(sttTable.getSTTtable().back().getNAnagenetic()==2);
        assert(sttTable.getSTTtable().back().getNColonisations()==1);
    }
    {  // STT update after merging archipelago does delete 'A' duplicates
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'A',false,4.0,4.0,{}),0);
        archi.addSpecies(Species(3.0,SpeciesID(1),SpeciesID(1),'A',true,4.0,4.0,{}),1);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateFullSTTtable(archi,3.0);
        assert(sttTable.getSTTtable().back().getNAnagenetic()==1);
        assert(sttTable.getSTTtable().back().getNColonisations()==1);
    }
    {  // STT update after merging archipelago does delete 'C' duplicates
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'C',false,4.0,4.0,{}),0);
        archi.addSpecies(Species(3.0,SpeciesID(1),SpeciesID(1),'C',true,4.0,4.0,{}),1);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateFullSTTtable(archi,3.0);
        assert(sttTable.getSTTtable().back().getNCladogenetic()==1);
        assert(sttTable.getSTTtable().back().getNColonisations()==1);
    }
    {  // STT update after merging archipelago does delete 'I' duplicates
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'I',false,4.0,3.9,{}),0);
        archi.addSpecies(Species(3.0,SpeciesID(1),SpeciesID(1),'I',true,4.0,3.9,{}),1);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateFullSTTtable(archi,3.0);
        assert(sttTable.getSTTtable().back().getNImmigrants()==1);
        assert(sttTable.getSTTtable().back().getNColonisations()==1);
    }
    {  // STT update after merging archipelago does delete 'I' duplicates
        // also when independent immigrations TODO: problem?
        int n_islands = 2;
        int islCarryingCap = 5;
        Archipelago archi = Archipelago(n_islands, islCarryingCap);
        archi.addSpecies(Species(3.9,SpeciesID(1),SpeciesID(1),'I',false,4.0,3.9,{}),0);
        archi.addSpecies(Species(3.0,SpeciesID(1),SpeciesID(1),'I',false,4.0,3.0,{}),1);
        STTtable sttTable = STTtable(1);
        assert(sttTable.size()==1);
        sttTable.updateFullSTTtable(archi,3.0);
        assert(sttTable.getSTTtable().back().getNImmigrants()==1);
        assert(sttTable.getSTTtable().back().getNColonisations()==1);
    }
}
