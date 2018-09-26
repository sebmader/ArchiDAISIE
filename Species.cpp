//
// Created by Sebastian Mader on 20.07.2018.
//

#include "Species.h"

Species::Species(const double birthTime, const int parentId,
        const int speciesId, const char status)
        : mBirthT{birthTime}, mParentID{parentId}, mSpeciesID{speciesId},
        mStatus{status}
{
    assert(birthTime >= 0.0);
    assert(parentId >= 0);
    assert(speciesId >= 0);
    assert(mStatus == 'I' || mStatus == 'A' || mStatus == 'C' || mStatus == 'M');
}

void Species::setBirth(const double time)
{
    mBirthT = time;
}

void Species::setStatus(char status)
{
    mStatus = status;
}


