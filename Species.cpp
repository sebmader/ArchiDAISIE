//
// Created by Sebastian Mader on 20.07.2018.
//

#include "Species.h"

Species::Species(const double birthTime, const SpeciesID parentId,
         const SpeciesID speciesId, const char status)
        : mBirthT{birthTime}, mParentID{parentId}, mSpeciesID{speciesId},
        mStatus{status}
{
    assert(birthTime >= 0.0);
    assert(mStatus == 'I' || mStatus == 'A' || mStatus == 'C' || mStatus == 'M'
            || mStatus == '0');
}

void Species::setBirth(const double time)
{
    mBirthT = time;
}

void Species::setStatus(char status)
{
    mStatus = status;
}

bool Species::isImmigrant() const noexcept
{
    return mStatus == 'I';
}

bool Species::isMigrant() const noexcept
{
    return mStatus == 'M';
}


