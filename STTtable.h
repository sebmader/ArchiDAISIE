//
// Created by Bastophiles on 08.11.2018.
//

#ifndef ARCHIDAISIE_STTTABLE_H
#define ARCHIDAISIE_STTTABLE_H

#include <vector>
#include <cassert>
#include "STT.h"
#include "Archipelago.h"
#include "Island.h"
#include "Species.h"

class STTtable {
public:
    explicit STTtable(const STT& firstRow = STT(0.0,0,0,0));
    void updateSTTtable(const Archipelago& archi, const double& time);

private:
    std::vector<STT> mSTTtable;
};

#endif //ARCHIDAISIE_STTTABLE_H
