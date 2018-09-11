//
// Created by Sebastian Mader on 05.09.2018.
//

#include "Island.h"
#include "Archipelago.h"
#include <vector>

using namespace std;
// ------ Island member functions ------

void Island::printIsland()
{
    for (auto &i : mvIsland)
        i.printSpec();
}

int Island::findPos(const int &ID) const    // find the position of certain species (input) in species vector
{                                                        // if not present: output = -1
    const int n_species = mvIsland.size();
    for (int i = 0; i < n_species; ++i)
        if (mvIsland[i].readSpID() == ID)
            return i;
    return -1;
}

int Island::createNewID()
{
    Archipelago::incrementMaxID();
    return Archipelago::returnMaxID();
}

double Island::calculateIslRates(
  const std::vector<double> &vIslPars, const int &iM, const int &iNumIsl,
        const double &dThisMLG)
{   // calculates the per-island rates of events, outputs them (for immidiate use)
    // and saves them within the island-class; input -> initial parameters
    // order of island parameter vector: gam_i, gam_m, lamb_cl, lamb_al, mu_l

    // calculate the per-island rates
    const int iN = static_cast<int>(mvIslSpecAlive.size());
    int iNImmi = 0;    // number of immigrant species -> only ones that can undergo anagenesis
    for (int g : mvIslSpecAlive)
        if (g <= 1000)
            ++iNImmi;
    const double dGam_I = max(0.0, (vIslPars[0] * iM * (1 - static_cast<double>(iN) / mIK)) / iNumIsl);  // immigration to specific island
    const double dGam_M = max(0.0, (vIslPars[1] * iN * dThisMLG) / iNumIsl*iNumIsl - iNumIsl); // migration from this island to all other islands
    // all other islands: '/iNumIsl', one-way: '/iNumIsl^2 - iNumIsl' -> ### CAUTION ### : really ???
    const double dLamb_Cl = max(0.0, vIslPars[2] * iN * (1 - (static_cast<double>(iN) / mIK))); // local cladogenesis
    const double dLamb_Al = max(0.0, vIslPars[3] * iNImmi);   // local anagenesis
    const double dMu_l = max(0.0, vIslPars[4] * iN);  // local extinction

    // save in rate vector
    mvLocalRates = { dGam_I, dGam_M, dLamb_Cl, dLamb_Al, dMu_l };

    double dSumIslRates = dGam_I + dGam_M + dLamb_Cl + dLamb_Al + dMu_l;

    return dSumIslRates;
}

double Island::extractSumIslRate() const noexcept
{
    double dSumRates = 0;
    for (double mvLocalRate : mvLocalRates)
        dSumRates += mvLocalRate;
    return dSumRates;
}

vector<int> Island::sampleLocalEvent(mt19937_64 prng, const int &iM)
{   // samples local event and species it happens to on THIS island
    // draw event
    const int iEvent = drawDisEvent(mvLocalRates, prng);
    // draw species
    int iSpecID;
    if (iEvent) {   // if not immigration (1-4) -> SpecID from extant island species
        iSpecID = mvIslSpecAlive[drawUniEvent(0, static_cast<int>(mvIslSpecAlive.size()), prng)];
    }
    else {  // if immigration (0) -> SpecID from mainland species pool
        iSpecID = drawUniEvent(1, iM, prng);
    }

    // return event and specID
    vector<int> vHappening = { iEvent, iSpecID };

    return vHappening;
}

// local updates:
void Island::immigrate(const int& iSpecID, double dTime)
{   // immigration from the mainland to THIS island
    // check if species is already present on island
    const int iPos = findPos(iSpecID);
    // if mainland sp -> BirthT = Time; else (if island sp) -> BirthT = old BirthT (because this function is used for migration as well !!
    Species newSpecies(dTime, iSpecID, iSpecID);
    if (iPos >= 0) {
        double newBirthT = mvIsland[iPos].readBirth();
        newSpecies = Species(newBirthT, iSpecID, iSpecID);    // immigrant
    }
    else {  // if present
        assert(iPos >= 0);
        assert(iPos < static_cast<int>(mvIsland.size()));
        if (!mvIsland[iPos].isExtant()) {   // if extinct
            mvIsland.push_back(newSpecies);
        }
        else  // if extant -> re-immigration ("re-setting the clock" (= BirthT))
            mvIsland[iPos] = newSpecies;
    }
    mvIslSpecAlive.push_back(iSpecID);
}

int Island::migrate(
  const int /* iSpecID */,
  vector<double> &vLogs,
  const double &dIniMigRate,
  mt19937_64 prng
)
{   // migration from THIS island to another; output: island of destination
    // draw island to which species migrates -> initial migration rate as parameter !!
    const int n_islands = vLogs.size();
    vector<double> vMigRates(n_islands);
    for (int k = 0; k < n_islands; ++k) {
        vMigRates[k] = max(0.0, (dIniMigRate * mvIslSpecAlive.size() * vLogs[k]) / n_islands*n_islands - n_islands);
    }
    const int iDestinationIsl = drawDisEvent(vMigRates, prng);

    return iDestinationIsl;
}

void Island::speciateClado(const int& iSpecID, double dTime)
{   // island species cladogenetically diverges
    // find species
    const int iPos = findPos(iSpecID);
    if (iPos == -1)
        throw logic_error("Drawn species is not present on island.. Something's wrong (Cladogenesis).\n");
    // 2 new species:
    Species spNew1 = Species(mvIsland[iPos].readBirth(), iSpecID, createNewID());
    Species spNew2 = Species(dTime, iSpecID, createNewID());
    // parent goes extinct
    mvIsland[iPos].goExtinct(dTime);
    mvIsland.push_back(spNew1);
    mvIsland.push_back(spNew2);
}

void Island::speciateAna(const int& iSpecID, double dTime)
{   // anagenetic speciation: only immigrants can undergo anagenesis !!!
    // find species
    const int iPos = findPos(iSpecID);
    if (iPos == -1)
        throw logic_error("Drawn species is not present on island.. Something's wrong (Anagenesis).\n");
    // new species
    Species spNew = Species(mvIsland[iPos].readBirth(), iSpecID, createNewID());
    // parent goes extinct & daugther gets added to island
    mvIsland[iPos].goExtinct(dTime);
    mvIsland.push_back(spNew);
}

void Island::goExtinct(const int& iSpecID, double dTime)
{   // species goes extinct
    // find species
    const int iPos = findPos(iSpecID);
    if (iPos == -1)
        throw logic_error("Drawn species is not present on island.. Something's wrong (Extinction).\n");
    // extinction
    mvIsland[iPos].goExtinct(dTime);
}

void Island::addIsland(const Island &islNew)
{   // adds another island to THIS (for aggregating archipelagos)
    // extract data frames from island that's to be added
    const vector<Species>& vSpecNew = islNew.returnIsland();
    const vector<int>& vSpecAliveNew = islNew.returnIslSpecAlive();

    // add species vector to THIS island
    if (!vSpecNew.empty()) {
        // intermediate vector of species
        vector<Species> vTempSpec = mvIsland;
        cerr << mvIsland.size() << '\n';
        cerr << vSpecNew.size() << '\n';
        vTempSpec.reserve(mvIsland.size() + vSpecNew.size());   // ### CAUTION ### : problem is the const type of ALL outputs and functions
        std::cerr << "HERE\n";                                                                        //  associated with .reserve()
        vTempSpec.insert(vTempSpec.end(), vSpecNew.begin(), vSpecNew.end());
        std::cerr << "PAST\n";
        // mvIsland.reserve(mvIsland.size() + vSpecNew.size());   // preallocate memory
        // mvIsland.insert(mvIsland.end(), vSpecNew.begin(), vSpecNew.end());
    }
    // add alive species to THIS island
    if (!vSpecAliveNew.empty()) {
        vector<int> vTempAlive = mvIslSpecAlive;
        vTempAlive.reserve(mvIslSpecAlive.size() + vSpecAliveNew.size());
        vTempAlive.insert(vTempAlive.end(), vSpecAliveNew.begin(), vSpecAliveNew.end());
        // mvIslSpecAlive.reserve(mvIslSpecAlive.size() + vSpecAliveNew.size());
        // mvIslSpecAlive.insert(mvIslSpecAlive.end(), vSpecAliveNew.begin(), vSpecAliveNew.end());
    }
}

const Species& Island::returnSpecies(const int pos) const
{
  assert(pos >= 0);
  assert(pos < static_cast<int>(mvIsland.size()));
  return mvIsland[pos];
}   // returns specific species from species vector
