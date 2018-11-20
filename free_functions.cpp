//
// Created by Sebastian Mader on 29.09.2018.
//

#include <fstream>
#include "free_functions.h"

using namespace std;

double getLogGrowth(const Island& island) noexcept
{
    return 1.0 - static_cast<double>(island.getNSpecies())
            / static_cast<double>(island.getCarryingCap());
}

double extractSumOfRates(const Island& island)
{
    vector<double> localRates = island.getLocalRates();
    double sumRates = 0;
    for (double rate : localRates)
        sumRates += rate;
    return sumRates;
}

Species findOldestSpecies(const vector<Species>& speciesVec)
{
    Species oldest = Species();
    for (auto& species : speciesVec) {
        if (species.getBirth() > oldest.getBirth())
            oldest = species;
    }
    return oldest;
}

STTtable mergeSTTtables(const std::vector<STTtable>& STTvec, const int& n_samples)
{
    double oldestAge = STTvec[0].getSTTtable()[0].getTime();
    for (auto& stt : STTvec) {
        assert(stt.getSTTtable()[0].getTime()==oldestAge);
    }
    STTtable fullSTT = STTtable((unsigned) n_samples+1, STT(oldestAge, 0, 0, 0, 0));

    // count species states per predefined time slice (island age / sample freq)
    for(size_t i = 1; i < fullSTT.size(); ++i) {
        double time = oldestAge - (oldestAge/n_samples) * i;
        assert(time >= 0.0);
        int nImmi = 0;
        int nAna = 0;
        int nClado = 0;
        int nColo = 0;
        for (auto& sttTab : STTvec) {
            for (int j = (int)sttTab.size()-1; j >= 0; --j) {
                if (sttTab[j].getTime() >= time) {
                    nImmi += sttTab[j].getNImmigrants();
                    nAna += sttTab[j].getNAnagenetic();
                    nClado += sttTab[j].getNCladogenetic();
                    nColo += sttTab[j].getNColonisations();
                    break;
                }
            }
        }
        fullSTT[i] = STT(time,nImmi,nAna,nClado,nColo);
    }
    return fullSTT;
}

int howManyLineages(const vector<Species>& species)
{
    vector<Species> lineages;
    for (auto& sp : species) {
        bool lineagePresent = false;
        for (auto& lineage : lineages) {
            if (lineage.isSister(sp))
                lineagePresent = true;
        }
        if(!lineagePresent)
            lineages.push_back(sp);
    }
    return static_cast<int>(lineages.size());
}

vector<SpeciesID> whichMainAncestors(const vector<Species>& species)
{
    vector<SpeciesID> mainlandAncestors;
    for (auto& sp : species) {
        bool lineagePresent = false;
        for (auto& mainAnc : mainlandAncestors) {
            if (mainAnc == sp.getParID())
                lineagePresent = true;
        }
        if(!lineagePresent)
            mainlandAncestors.push_back(sp.getParID());
    }
    return mainlandAncestors;
}

void outputBranching(const Island& fullIsland, ofstream& ofs)
{  // output: datatable(clade name, status, missing species, branching times),
    // island age,
    // mainland species
    const vector<Species>& species = fullIsland.getSpecies();
    if(!ofs.is_open()) {
        throw runtime_error("unable to open file.");
    }
    if (species.empty()) {
        // TODO: what if archipelago is empty? -> return output with empty data table
        ofs << "" << ',' << "" << ',' << "" << ',' << "" << ',' << '\n';
    }
    else {
        // data table
          // how many clades -> size of datatable
          // collect species of each clade
        vector<SpeciesID> mainAncestors = whichMainAncestors(species);
        const int n_mainColoniser = static_cast<int>(mainAncestors.size());
        vector<vector<Species> > clades(n_mainColoniser, vector<Species>());
        for (int i = 0; i < n_mainColoniser; ++i) {
            for (auto& sp : species) {
                if (sp.getParID() == mainAncestors[i]) {
                    clades[i].push_back(sp);
                }
            }
        }
        // extract status -> endemic, nonendemic or both?
        vector<string> status(mainAncestors.size());
        for (int i = 0; i < n_mainColoniser; ++i) {
            bool isEndemic = false;
            bool isNonEndemic = false;
            for (auto& sp : clades[i]) {
                if (sp.getStatus() == 'I')
                    isNonEndemic = true;
                else if (sp.getStatus() == 'C' || sp.getStatus() == 'A')
                    isEndemic = true;
                else
                    throw logic_error("unknown status");
            }
            if (isEndemic && isNonEndemic)
                status[i] = "Endemic&Non_Endemic";
            else if (isEndemic && !isNonEndemic)
                status[i] = "Endemic";
            else if (!isEndemic && isNonEndemic)
                status[i] = "Non_endemic";
            else
                throw logic_error("unknown status");
        }
        // extract birth times -> make string
        vector<string> branchingTimes((unsigned)n_mainColoniser);
        for (int i = 0; i < n_mainColoniser; ++i) {
            for (auto& sp : clades[i]) {
                branchingTimes[i] += to_string(sp.getBirth()) + ',';
            }
            branchingTimes[i].pop_back(); // no comma at the end
        }
        // names == mainland ancestor
        // output to file
        for (int k = 0; k < n_mainColoniser; ++k) {
            ofs << mainAncestors[k].getSpeciesID() << ',' << '\"'+status[k]+'\"' << ','
                << 0 << ',' << '\"'+branchingTimes[k]+'\"' << '\n';
        }
    }
}


void outputBranching(const Island& fullIsland, ostream& os)
{  // output: datatable(clade name, status, missing species, branching times),
    // island age,
    // mainland species
    const vector<Species>& species = fullIsland.getSpecies();
    if (species.empty()) {
        // TODO: what if archipelago is empty? -> return output with empty data table
        os << "" << ',' << "" << ',' << "" << ',' << "" << ',' << '\n';
    }
    else {
        // data table
        // how many clades -> size of datatable
        // collect species of each clade
        vector<SpeciesID> mainAncestors = whichMainAncestors(species);
        const int n_mainColoniser = static_cast<int>(mainAncestors.size());
        vector<vector<Species> > clades(n_mainColoniser, vector<Species>());
        for (int i = 0; i < n_mainColoniser; ++i) {
            for (auto& sp : species) {
                if (sp.getParID() == mainAncestors[i]) {
                    clades[i].push_back(sp);
                }
            }
        }
        // extract status -> endemic, nonendemic or both?
        vector<string> status(mainAncestors.size());
        for (int i = 0; i < n_mainColoniser; ++i) {
            bool isEndemic = false;
            bool isNonEndemic = false;
            for (auto& sp : clades[i]) {
                if (sp.getStatus() == 'I')
                    isNonEndemic = true;
                else if (sp.getStatus() == 'C' || sp.getStatus() == 'A')
                    isEndemic = true;
                else
                    throw logic_error("unknown status");
            }
            if (isEndemic && isNonEndemic)
                status[i] = "Endemic&Non_Endemic";
            else if (isEndemic && !isNonEndemic)
                status[i] = "Endemic";
            else if (!isEndemic && isNonEndemic)
                status[i] = "Non_endemic";
            else
                throw logic_error("unknown status");
        }
        // extract birth times -> make string
        vector<string> branchingTimes((unsigned)n_mainColoniser);
        for (int i = 0; i < n_mainColoniser; ++i) {
            for (auto& sp : clades[i]) {
                branchingTimes[i] += to_string(sp.getBirth()) + ',';
            }
            branchingTimes[i].pop_back(); // no comma at the end
        }
        // names == mainland ancestor
        // output to file
        for (int k = 0; k < n_mainColoniser; ++k) {
            os << mainAncestors[k].getSpeciesID() << ',' << '\"'+status[k]+'\"' << ','
                << 0 << ',' << '\"'+branchingTimes[k]+'\"' << '\n';
        }
    }
}