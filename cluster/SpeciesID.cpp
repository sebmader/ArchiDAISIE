//
// Created by Bastophiles on 12.09.2018.
//

#include <cassert>
#include "SpeciesID.h"

SpeciesID::SpeciesID(const int id) : speciesID{id} {
    assert(id >= 0);
}

int SpeciesID::getSpeciesID() const noexcept
{
    return speciesID;
}

SpeciesID SpeciesID::createNewSpeciesID()
{
    incrementSpeciesID();
    return static_cast<SpeciesID>(speciesID);
}

bool SpeciesID::operator==(const SpeciesID& ID1) const
{
    return speciesID==ID1.speciesID;
}

bool SpeciesID::operator!=(const SpeciesID& ID1) const
{
    return !(ID1==*this);
}


bool SpeciesID::operator<(const SpeciesID& rhs) const
{
    return speciesID<rhs.speciesID;
}

bool SpeciesID::operator>(const SpeciesID& rhs) const
{
    return rhs<*this;
}

bool SpeciesID::operator<=(const SpeciesID& rhs) const
{
    return !(rhs<*this);
}

bool SpeciesID::operator>=(const SpeciesID& rhs) const
{
    return !(*this<rhs);
}

