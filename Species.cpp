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
        const double colonisationT,
        const vector<char> cladoStates)
        : mBirthT{birthTime}, mMainParentID{parentId}, mSpeciesID{speciesId},
        mStatus{status}, mMigrated{migrated}, mAncestralBT{ancestralBT},
        mColonisationT{colonisationT}, mCladoStates{cladoStates}
{
    assert(birthTime >= 0.0);
    assert(status == 'I' || status == 'A' || status == 'C' || status == '0');
    assert(ancestralBT >= 0.0);
    assert(colonisationT >= 0.0);
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

bool Species::isSister(const Species& potentialSis) const
{
    if(*this == potentialSis)
        return false;
    return mMainParentID == potentialSis.mMainParentID &&
                mColonisationT == potentialSis.mColonisationT;
}

bool Species::isMostRecentSis(const Species& potentialSis) const
{
    vector<char> formerSpeciations = mCladoStates;
    formerSpeciations.pop_back();
    vector<char> sisFormerSpeciations = potentialSis.getCladoStates();
    while(sisFormerSpeciations.size() > formerSpeciations.size())
        sisFormerSpeciations.pop_back();

    if (!this->isSister(potentialSis) || formerSpeciations != sisFormerSpeciations)
        return false;

    const int posLastSpeciation = static_cast<int>(mCladoStates.size() - 1);
    char sisLastCladoState = 'a';
    if (mCladoStates.back() == 'a')
        sisLastCladoState = 'b';
    return potentialSis.mCladoStates[posLastSpeciation] == sisLastCladoState;
}

void Species::printSpec()
{
    std::cout << mBirthT << '\t' << mMainParentID.getSpeciesID() << '\t'
              << mSpeciesID.getSpeciesID() << '\t' << mStatus << '\t' << mMigrated << '\t'
              << mAncestralBT << '\t' << mColonisationT << '\t';
    for(auto& cladoState : mCladoStates)
        cout << cladoState;
    cout <<  '\n';
}

bool Species::operator==(const Species& rhs) const
{
    return mBirthT==rhs.mBirthT &&
            mMainParentID==rhs.mMainParentID &&
            mSpeciesID==rhs.mSpeciesID &&
            mStatus==rhs.mStatus &&
            mMigrated==rhs.mMigrated &&
            mAncestralBT==rhs.mAncestralBT &&
            mColonisationT==rhs.mColonisationT &&
            mCladoStates==rhs.mCladoStates;
}

bool Species::operator!=(const Species& rhs) const
{
    return !(rhs==*this);
}

