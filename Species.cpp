//
// Created by Sebastian Mader on 20.07.2018.
//

#include "Species.h"

Species::Species(const double birthT, const int parentId, const int speciesId)
        : dBirthT{birthT}, iParentId{parentId}, iSpecId{speciesId}, dExtinctT{-1.0}
{
    assert(dBirthT >= 0.0);
    assert(iParentId >= 0);
    assert(iSpecId >= 0);
    assert(isExtant());
}

bool Species::isExtant() const noexcept
{
    return dExtinctT == -1.0;
}


