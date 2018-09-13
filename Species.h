//
// Created by Sebastian Mader on 20.07.2018.
//

#ifndef ARCHIDAISIE_SPECIES_H
#define ARCHIDAISIE_SPECIES_H

#include <iostream>
#include <cassert>


class Species {

public:
    Species(double birthT, int parentId, int speciesId);

    void goExtinct(const double time) {dExtinctT = time;}
    bool isExtant() const noexcept;
    double readBirth() const noexcept {return dBirthT;}
    int readParID() const noexcept {return iParentId;}
    int readSpID() const noexcept {return iSpecId;}
    const double& readExtinct() const noexcept {return dExtinctT;}

    void printSpec() {
        std::cout << dBirthT << '\t' << iParentId << '\t'
                  << iSpecId << '\t' << dExtinctT << '\n';
    }

private:
    double dBirthT;     //Should be const one day
    int iParentId;      //Should be const one day
    int iSpecId;        //Should be const one day
    double dExtinctT;   //Will be -1 when extant
};



#endif //ARCHIDAISIE_SPECIES_H
