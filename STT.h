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
            int nCladogenetic = 0,
            int nColonisations = 0);

    double getTime() const noexcept;
    int getNImmigrants() const noexcept;
    int getNAnagenetic() const noexcept;
    int getNCladogenetic() const noexcept;
    int getNColonisations() const noexcept;

private:
    double mTime;
    int mNImmigrants;
    int mNAnagenetic;
    int mNCladogenetic;
    int mNColonisations;
};

#endif //ARCHIDAISIE_STT_H
