//
// Created by Bastophiles on 12.09.2018.
//

#ifndef ARCHIDAISIE_SPECIESID_H
#define ARCHIDAISIE_SPECIESID_H

class SpeciesID {
public:
    explicit SpeciesID(int = 0);

    bool operator==(const SpeciesID& ID1) const;
    bool operator!=(const SpeciesID& ID1) const;
    bool operator<(const SpeciesID& rhs) const;
    bool operator>(const SpeciesID& rhs) const;
    bool operator<=(const SpeciesID& rhs) const;
    bool operator>=(const SpeciesID& rhs) const;


    void incrementSpeciesID() noexcept { ++speciesID; }
    int getSpeciesID() const noexcept;

    SpeciesID createNewSpeciesID();

private:
    int speciesID;
};

#endif //ARCHIDAISIE_SPECIESID_H
