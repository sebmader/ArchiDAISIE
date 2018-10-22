//
// Created by Sebastian Mader on 20.07.2018.
//

#include "Species.h"
using namespace std;

Species::Species(const double birthTime,
        const SpeciesID parentId,
        const SpeciesID speciesId,
        const char status,
        const bool migrated,
        const double ancestralBT,
        const vector<char> cladoStates)
        : mBirthT{birthTime}, mMainParentID{parentId}, mSpeciesID{speciesId},
        mStatus{status}, mMigrated{migrated}, mAncestralBT{ancestralBT},
        mCladoStates{cladoStates}
{
    assert(birthTime >= 0.0);
    assert(status == 'I' || status == 'A' || status == 'C' || status == '0');
    assert(ancestralBT >= 0.0);
}

void Species::setBirth(const double time)
{
    assert(time >= 0.0);
    mBirthT = time;
}

void Species::setStatus(char status)
{
    assert(status == 'I' || status == 'A' || status == 'C' || status == 'M'
            || status == '0');
    mStatus = status;
}

void Species::setAncestralBT(const double time)
{
    assert(time >= 0.0);
    mAncestralBT = time;
}

bool Species::isImmigrant() const noexcept
{
    return mStatus == 'I';
}

bool Species::hasMigrated() const noexcept
{
    return mMigrated;
}

bool Species::isCladogenetic() const noexcept
{
    return mStatus == 'C';
}


