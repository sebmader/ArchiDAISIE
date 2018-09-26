//
// Created by Sebastian Mader on 20.07.2018.
//

#include "Species.h"

Species::Species(const double birthTime, const int parentId, const int speciesId)
        : birthT{birthTime}, parentID{parentId}, speciesID{speciesId}, extinctT{-1.0}
{
    assert(birthTime >= 0.0);
    assert(parentId >= 0);
    assert(speciesId >= 0);
    assert(isExtant());
}

bool Species::isExtant() const noexcept
{
    return extinctT == -1.0;
}

void Species::setBirth(const double time)
{
    birthT = time;
}


