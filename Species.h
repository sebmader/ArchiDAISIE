//
// Created by Sebastian Mader on 20.07.2018.
//

#ifndef ARCHIDAISIE_SPECIES_H
#define ARCHIDAISIE_SPECIES_H

#include <iostream>
#include <cassert>
#include <vector>
#include "SpeciesID.h"


class Species {

public:
    explicit Species(double = 0.0, SpeciesID = SpeciesID(), SpeciesID = SpeciesID(),
            char = '0', double = 0.0, std::vector<char> = {});

    double getBirth() const noexcept { return mBirthT; }
    SpeciesID getSpecID() const noexcept { return mSpeciesID; }
    SpeciesID getParID() const noexcept { return mMainParentID; }
    char getStatus() const noexcept { return mStatus; }
    double getAncestralBT() const noexcept { return mAncestralBT; }
    std::vector<char> getCladoStates() const noexcept { return mCladoStates; }
    void setBirth(double);
    void setStatus(char);
    void setAncestralBT(double);
    bool isImmigrant() const noexcept;
    bool isMigrant() const noexcept;
    bool isCladogenetic() const noexcept;



    void printSpec() {
        std::cout << mBirthT << '\t' << mMainParentID.getSpeciesID() << '\t'
                  << mSpeciesID.getSpeciesID() << '\t' << mStatus
                  << '\t' << mAncestralBT << '\n';
    }

private:
    double mBirthT;     // time of birth/colonisation of species
    SpeciesID mMainParentID;      // ID of parent species TODO: test Mainland parent !!
    SpeciesID mSpeciesID;     // ID of species
    char mStatus;       // immigrant, anagenesis, cladogenesis, migrant, default
                            // ('I','A','C','M','0')
    double mAncestralBT;     // birth time of ancestor = colonisation/cladogenesis (if 2nd daughter)
                            // saved before migration (gives species a new birth time)
                            // -> need to save this for preparing data for DAISIE_ML
                            // is equal to mBirthT if hasn't migrated
    std::vector<char> mCladoStates;     // vector of states after speciation: whether it is
                                            // daughter 'A' or daughter 'B'
};


#endif //ARCHIDAISIE_SPECIES_H
