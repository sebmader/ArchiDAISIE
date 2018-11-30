//
// Created by Sebastian Mader on 08.11.2018.
//

#include "STT.h"

STT::STT(const double time, const int& nImmigrants, const int& nAnagenetic, const int& nCladogenetic,
        const int& nColonisations)
        : mTime(time), mNImmigrants(nImmigrants), mNAnagenetic(nAnagenetic),
        mNCladogenetic(nCladogenetic), mNColonisations(nColonisations)
{
    assert(time >= 0.0);
    assert(mTime == time);
    assert(nImmigrants >= 0);
    assert(nAnagenetic >= 0);
    assert(nCladogenetic >= 0);
    assert(nColonisations >= 0);
}

double STT::getTime() const noexcept
{
    return mTime;
}

int STT::getNImmigrants() const noexcept
{
    return mNImmigrants;
}

int STT::getNAnagenetic() const noexcept
{
    return mNAnagenetic;
}

int STT::getNCladogenetic() const noexcept
{
    return mNCladogenetic;
}

int STT::getNColonisations() const noexcept
{
    return mNColonisations;
}

bool STT::operator==(const STT& rhs) const
{
    return mTime==rhs.mTime &&
            mNImmigrants==rhs.mNImmigrants &&
            mNAnagenetic==rhs.mNAnagenetic &&
            mNCladogenetic==rhs.mNCladogenetic &&
            mNColonisations==rhs.mNColonisations;
}

bool STT::operator!=(const STT& rhs) const
{
    return !(rhs==*this);
}
