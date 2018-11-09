//
// Created by Sebastian Mader on 08.11.2018.
//

#include "STT.h"

STT::STT(const double& time, int nImmigrants, int nAnagenetic, int nCladogenetic, int nColonisations)
        :mTime(time), mNImmigrants(nImmigrants), mNAnagenetic(nAnagenetic),
        mNCladogenetic(nCladogenetic)
{
    assert(time >= 0.0);
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
