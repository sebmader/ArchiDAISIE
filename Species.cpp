//
// Created by Sebastian Mader on 20.07.2018.
//

#include "Species.h"
using namespace std;

Species::Species(const double birthTime, const SpeciesID parentId,
         const SpeciesID speciesId, const char status, const double cladeBirth,
         const vector<char> cladoStates)
        : mBirthT{birthTime}, mParentID{parentId}, mSpeciesID{speciesId},
        mStatus{status}, mCladeBirthT{cladeBirth}, mCladoStates{cladoStates}
{
    assert(birthTime >= 0.0);
    assert(status == 'I' || status == 'A' || status == 'C' || status == 'M'
            || status == '0');
    assert(cladeBirth >= 0.0);
}

void Species::setBirth(const double time)
{
    mBirthT = time;
}

void Species::setStatus(char status)
{
    mStatus = status;
}

void Species::setCladeBirth(const double time)
{
    mCladeBirthT = time;
}

bool Species::isImmigrant() const noexcept
{
    return mStatus == 'I';
}

bool Species::isMigrant() const noexcept
{
    return mStatus == 'M';
}

bool Species::isCladogenetic() const noexcept
{
    return mStatus == 'C';
}


