//
// Created by Sebastian Mader on 08.11.2018.
//

#ifndef ARCHIDAISIE_STTTABLE_H
#define ARCHIDAISIE_STTTABLE_H

#include <vector>
#include <cassert>
#include <ostream>
#include "STT.h"
#include "Archipelago.h"
#include "Species.h"

class STTtable {
public:
    explicit STTtable(const unsigned int& size = 0,
            const STT& rowContent = STT());

    const std::vector<STT>& getSTTtable() const noexcept;
    unsigned long size() const;

    void push_back(const STT& newRow);
    void updateSTTtable(const Archipelago& archi, const double& time);

    friend std::ostream& operator<<(std::ostream& os, const STTtable& table);

private:
    std::vector<STT> mSTTtable;
};

#endif //ARCHIDAISIE_STTTABLE_H
