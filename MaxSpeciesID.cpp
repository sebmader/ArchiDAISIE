//
// Created by Bastophiles on 12.09.2018.
//

#include <cassert>
#include "MaxSpeciesID.h"

MaxSpeciesID::MaxSpeciesID(const int maxID) : maxSpeciesID{maxID} {
    assert(maxID >= 0);
}

int MaxSpeciesID::createNewSpeciesID()
{
    incrementMaxSpeciesID();
    return getMaxSpeciesID();
}

