//
// Created by Sebastian Mader on 20.07.2018.
//

#ifndef ARCHIDAISIE_SPECIES_H
#define ARCHIDAISIE_SPECIES_H


class Species {

public:
    Species(const double birthT, const int parentId, const int speciesId);

    void goExtinct(const double &time) {dExtinctT = time;}
    bool isExtant() const noexcept {return dExtinctT == -1.0;}
    const double& readBirth() const noexcept {return dBirthT;}
    const int& readParID() const noexcept {return iParentId;}
    const int& readSpID() const noexcept {return iSpecId;}
    const double& readExtinct() const noexcept {return dExtinctT;}

    void printSpec() {
        std::cout << dBirthT << '\t' << iParentId << '\t' << iSpecId << '\t' << dExtinctT << std::endl;
    }

private:
    double dBirthT;     //Should be const one day
    int iParentId;      //Should be const one day
    int iSpecId;        //Should be const one day
    double dExtinctT;   //Will be -1 when extant
};

Species::Species(const double birthT, const int parentId, const int speciesId)
        : dBirthT{birthT}, iParentId{parentId}, iSpecId{speciesId}, dExtinctT{-1.0}
{
    assert(dBirthT >= 0.0);
    assert(iParentId >= 0);
    assert(iSpecId >= 0);
    assert(isExtant());
}

#endif //ARCHIDAISIE_SPECIES_H
