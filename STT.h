//
// Created by Sebastian Mader on 08.11.2018.
//

#ifndef ARCHIDAISIE_STT_H
#define ARCHIDAISIE_STT_H

#include <cassert>

class STT {
public:
    explicit STT(const double& time = 0.0,
            int nImmigrants = 0,
            int nAnagenetic = 0,
            int nCladogenetic = 0);

    double getTime() const;

    int getNImmigrants() const;

    int getNAnagenetic() const;

    int getNCladogenetic() const;

private:
    const double mTime;
    const int mNImmigrants;
    const int mNAnagenetic;
    const int mNCladogenetic;
};

#endif //ARCHIDAISIE_STT_H
