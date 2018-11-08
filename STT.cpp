//
// Created by Sebastian Mader on 08.11.2018.
//

#include "STT.h"

STT::STT(const double& time, const int nImmigrants, const int nAnagenetic, const int nCladogenetic)
        :mTime(time), mNImmigrants(nImmigrants), mNAnagenetic(nAnagenetic),
        mNCladogenetic(nCladogenetic)
{
    assert(time >= 0.0);
    assert(nImmigrants >= 0);
    assert(nAnagenetic >= 0);
    assert(nCladogenetic >= 0);
}

double STT::getTime() const
{
    return mTime;
}

int STT::getNImmigrants() const
{
    return mNImmigrants;
}

int STT::getNAnagenetic() const
{
    return mNAnagenetic;
}

int STT::getNCladogenetic() const
{
    return mNCladogenetic;
}
