//
// Created by Sebastian Mader on 20.07.2018.
//

#ifndef ARCHIDAISIE_SPECIES_H
#define ARCHIDAISIE_SPECIES_H

#include <iostream>
#include <cassert>
#include "SpeciesID.h"


class Species {

public:
    explicit Species(double = 0, SpeciesID = SpeciesID(), SpeciesID = SpeciesID(), char = 'I');

    double readBirth() const noexcept { return mBirthT; }
    SpeciesID readSpID() const noexcept { return mSpeciesID; }
    char readStat() const noexcept { return mStatus; }
    void setBirth(double);
    void setStatus(char);
    bool isImmigrant() const noexcept;
    bool isMigrant() const noexcept;

    void printSpec() {
        std::cout << mBirthT << '\t' << mParentID.getMaxSpeciesID() << '\t'
                  << mSpeciesID.getMaxSpeciesID() << '\t' << mStatus << '\n';
    }

private:
    double mBirthT;     //Should be const one day
    SpeciesID mParentID;      //Should be const one day
    SpeciesID mSpeciesID;     //Should be const one day
    char mStatus;       // immigrant, anagenesis, cladogenesis, migrant ('I','A','C','M')
};


#endif //ARCHIDAISIE_SPECIES_H
