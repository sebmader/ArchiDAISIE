//
// Created by Sebastian Mader on 08.11.2018.
//

#ifndef ARCHIDAISIE_STTTABLE_H
#define ARCHIDAISIE_STTTABLE_H

#include <vector>
#include <cassert>
#include <ostream>
#include "STT.h"
#include "SpeciesID.h"
#include "Species.h"
#include "Island.h"
#include "Archipelago.h"


class STTtable {
public:
    explicit STTtable(const int size = 0);

    const std::vector<STT>& getSTTtable() const noexcept;

    unsigned long size();
    void updateSTTtable(const Archipelago& archi, const double& time);

    friend std::ostream& operator<<(std::ostream& os, const STTtable& table);

private:
    std::vector<STT> mSTTtable;
};

#endif //ARCHIDAISIE_STTTABLE_H
