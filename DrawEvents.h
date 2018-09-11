//
// Created by Bastophiles on 06.09.2018.
//

#ifndef ARCHIDAISIE_DRAWEVENTS_H
#define ARCHIDAISIE_DRAWEVENTS_H

#include <vector>
#include <iostream>
#include <random>

int drawDisEvent(const std::vector<double> &, std::mt19937_64&);

int drawUniEvent(const int &, const int &, std::mt19937_64&);

#endif //ARCHIDAISIE_DRAWEVENTS_H
