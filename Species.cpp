//
// Created by Sebastian Mader on 20.07.2018.
//

#include "Species.h"

Species::Species(const double birthT, const int parentId, const int speciesId)
        : birthT{birthT}, parentID{parentId}, speciesID{speciesId}, extinctT{-1.0}
{
    assert(birthT >= 0.0);
    assert(parentId >= 0);
    assert(speciesID >= 0);
    assert(isExtant());
}

bool Species::isExtant() const noexcept
{
    return extinctT == -1.0;
}


