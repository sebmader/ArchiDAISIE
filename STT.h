//
// Created by Sebastian Mader on 08.11.2018.
//

#ifndef ARCHIDAISIE_STT_H
#define ARCHIDAISIE_STT_H

#include <cassert>

class STT {
public:
    STT(const double& time, int nImmigrants, int nAnagenetic, int nCladogenetic);
private:
    const double mTime;
    const int mNImmigrants;
    const int mNAnagenetic;
    const int mNCladogenetic;
};

#endif //ARCHIDAISIE_STT_H
