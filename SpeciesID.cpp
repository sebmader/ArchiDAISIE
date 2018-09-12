//
// Created by Bastophiles on 12.09.2018.
//

#include <cassert>
#include "SpeciesID.h"

SpeciesID::SpeciesID(const int maxID) : maxSpeciesID{maxID} {
    assert(maxID >= 0);
}

int SpeciesID::createNewSpeciesID()
{
    incrementMaxSpeciesID();
    return getMaxSpeciesID();
}

