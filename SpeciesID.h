//
// Created by Bastophiles on 12.09.2018.
//

#ifndef ARCHIDAISIE_SPECIESID_H
#define ARCHIDAISIE_SPECIESID_H

class SpeciesID {
public:
    explicit SpeciesID(int = 0);

    void incrementMaxSpeciesID() noexcept { ++maxSpeciesID; }
    int getMaxSpeciesID() const noexcept { return maxSpeciesID; }

    int createNewSpeciesID();

private:
    int maxSpeciesID;
};

#endif //ARCHIDAISIE_SPECIESID_H
