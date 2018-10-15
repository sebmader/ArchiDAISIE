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
    explicit Species(double = 0.0, SpeciesID = SpeciesID(), SpeciesID = SpeciesID(),
            char = '0', double = 0.0);

    double getBirth() const noexcept { return mBirthT; }
    SpeciesID getSpecID() const noexcept { return mSpeciesID; }
    SpeciesID getParID() const noexcept { return mParentID; }
    char getStatus() const noexcept { return mStatus; }
    double getCladeBirthT() const noexcept { return mCladeBirthT; }
    void setBirth(double);
    void setStatus(char);
    void setCladeBirth(double);
    bool isImmigrant() const noexcept;
    bool isMigrant() const noexcept;



    void printSpec() {
        std::cout << mBirthT << '\t' << mParentID.getSpeciesID() << '\t'
                  << mSpeciesID.getSpeciesID() << '\t' << mStatus
                  << '\t' << mCladeBirthT << '\n';
    }

private:
    double mBirthT;     // Should be const one day
    SpeciesID mParentID;      // Should be const one day
    SpeciesID mSpeciesID;     // Should be const one day
    char mStatus;       // immigrant, anagenesis, cladogenesis, migrant, default
                            // ('I','A','C','M','0')
    double mCladeBirthT;     // birth time of clade = colonisation/cladogenesis (if 2nd daughter)
                            // saved before migration (gives species a new birth time)
                            // -> need to save this for preparing data for DAISIE_ML
                            // is equal to mBirthT if hasn't migrated yet
};


#endif //ARCHIDAISIE_SPECIES_H
