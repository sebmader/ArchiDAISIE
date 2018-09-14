//
// Created by Sebastian Mader on 20.07.2018.
//

#ifndef ARCHIDAISIE_SPECIES_H
#define ARCHIDAISIE_SPECIES_H

#include <iostream>
#include <cassert>


class Species {

public:
    explicit Species(double = 0, int = 0, int = 0);

    void goExtinct(const double time) {extinctT = time;}
    bool isExtant() const noexcept;
    double readBirth() const noexcept {return birthT;}
    int readParID() const noexcept {return parentID;}
    int readSpID() const noexcept {return speciesID;}
    const double& readExtinct() const noexcept {return extinctT;}

    void printSpec() {
        std::cout << birthT << '\t' << parentID << '\t'
                  << speciesID << '\t' << extinctT << '\n';
    }

private:
    double birthT;     //Should be const one day
    int parentID;      //Should be const one day
    int speciesID;        //Should be const one day
    double extinctT;   //Will be -1 when extant
};



#endif //ARCHIDAISIE_SPECIES_H
