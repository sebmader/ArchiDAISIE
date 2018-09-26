//
// Created by Sebastian Mader on 20.07.2018.
//

#ifndef ARCHIDAISIE_SPECIES_H
#define ARCHIDAISIE_SPECIES_H

#include <iostream>
#include <cassert>


class Species {

public:
    explicit Species(double = 0, int = 0, int = 0, char = 'I');

    double readBirth() const noexcept { return mBirthT; }
    int readSpID() const noexcept { return mSpeciesID; }
    char readStat() const noexcept { return mStatus; }
    void setBirth(double);
    void setStatus(char);

    void printSpec() {
        std::cout << mBirthT << '\t' << mParentID << '\t'
                  << mSpeciesID << '\t' << mStatus << '\n';
    }

private:
    double mBirthT;     //Should be const one day
    int mParentID;      //Should be const one day
    int mSpeciesID;     //Should be const one day
    char mStatus;       // immigrant, anagenesis, cladogenesis, migrant ('I','A','C','M')
};


#endif //ARCHIDAISIE_SPECIES_H
