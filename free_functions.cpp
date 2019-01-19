//
// Created by Sebastian Mader on 29.09.2018.
//

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
    double sumRates = 0.0;
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

vector<Species> whichLineageAncestors(const vector<Species>& species)
{  // outputs the reconstructed mainland ancestor (BirthTime, ID) of each lineage
    // could reconstruct whole species, but then it wouldn't find itself anymore if extant
    vector<Species> lineageExample;
    for (auto& sp : species) {
        bool lineagePresent = false;
        for (auto& lineage : lineageExample) {
            if (lineage.isSister(sp))
                lineagePresent = true;
        }
        if(!lineagePresent)
            lineageExample.push_back(sp);
    }
    vector<Species> lineageAncestors(lineageExample.size());
    for (size_t i = 0; i < lineageExample.size(); ++i) {
        Species mainAncestor = Species(lineageExample[i].getColonisationT(),
                lineageExample[i].getParID(),
                lineageExample[i].getParID(),
                '0', false, 0.0,
                lineageExample[i].getColonisationT());
        lineageAncestors[i] = mainAncestor;
    }
    return lineageAncestors;
}

vector<Species> whichMainAncestors(const vector<Species>& species)
{  // outputs the reconstructed mainland ancestor (BirthTime, ID) of each lineage
    // could reconstruct whole species, but then it wouldn't find itself anymore if extant
    vector<Species> examplePerMainSpec;
    for (auto& sp : species) {
        bool lineagePresent = false;
        for (auto& lineage : examplePerMainSpec) {
            if (lineage.getParID() == sp.getParID()) {
                lineagePresent = true;
                if (lineage.getColonisationT() < sp.getColonisationT())
                    lineage = sp;
            }
        }
        if(!lineagePresent)
            examplePerMainSpec.push_back(sp);
    }
    vector<Species> mainlandAncestors(examplePerMainSpec.size());
    for (size_t i = 0; i < examplePerMainSpec.size(); ++i) {
        Species mainAncestor = Species(examplePerMainSpec[i].getColonisationT(),
                examplePerMainSpec[i].getParID(),
                examplePerMainSpec[i].getParID(),
                '0', false, 0.0,
                examplePerMainSpec[i].getColonisationT());
        mainlandAncestors[i] = mainAncestor;
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
    string column_names = "Clade_name,Status,Missing_species,Branching_times";
    ofs << column_names << '\n';
    if (species.empty()) {
        ofs << "" << ',' << "" << ',' << "" << ',' << "" << ',' << '\n';
    }
    else {
        // data table
          // how many clades -> size of datatable
          // collect species of each clade
        vector<Species> mainAncestors = whichMainAncestors(species);
        const unsigned n_mainColoniser = static_cast<int>(mainAncestors.size());
        vector<vector<Species> > missingSpecies(n_mainColoniser, vector<Species>());
        vector<vector<Species> > clades(n_mainColoniser, vector<Species>());
        for (size_t i = 0; i < n_mainColoniser; ++i) {
            for (auto& sp : species) {
                if (sp.isSister(mainAncestors[i])) {
                    clades[i].push_back(sp);
                }
                // re-immigrants = missing species
                else if (sp.getSpecID() == mainAncestors[i].getSpecID()) {
                    missingSpecies[i].push_back(sp);
                }
            }
        }
        // extract status -> endemic, nonendemic or both?
        vector<string> status(mainAncestors.size());
        for (size_t i = 0; i < n_mainColoniser; ++i) {
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
            for (auto& missSp : missingSpecies[i]) {
                if (missSp.getStatus() == 'I')
                    isNonEndemic = true;
                else if (missSp.getStatus() == 'C' || missSp.getStatus() == 'A')
                    isEndemic = true;
                else
                    throw logic_error("unknown status");
            }
            if (isEndemic && isNonEndemic)
                status[i] = "Endemic&Non_endemic";
            else if (isEndemic && !isNonEndemic)
                status[i] = "Endemic";
            else if (!isEndemic && isNonEndemic)
                status[i] = "Non_endemic";
            else
                throw logic_error("unknown status");
        }
        // extract birth times -> make string
        vector<string> branchingTimes(n_mainColoniser);
        for (size_t i = 0; i < n_mainColoniser; ++i) {
            for (auto& sp : clades[i]) {
                branchingTimes[i] += to_string(sp.getBirth()) + ',';
            }
            branchingTimes[i].pop_back(); // no comma at the end
        }
        // names == mainland ancestor
        // output to file
        for (size_t k = 0; k < n_mainColoniser; ++k) {
            ofs << mainAncestors[k].getSpecID().getSpeciesID() << ',' << '\"'+status[k]+'\"' << ','
                << missingSpecies[k].size() << ',' << '\"'+branchingTimes[k]+'\"' << '\n';
        }
    }
}

vector<string> readParameterRowCSV(ifstream& inputFile, const int& lineWanted)
{
    if(!inputFile.is_open())
        throw runtime_error("unable to open input file.\n");

    int lineNum = 0;
    string row;
    vector<string> splittedRowCSV;
    while (getline(inputFile,row)) {
        ++lineNum;
        if (lineNum == lineWanted) {
            string cell;
            istringstream buffer(row);
            while (getline(buffer, cell, ',')) {
                splittedRowCSV.push_back(cell);
            }
        }
    }
    return splittedRowCSV;
}

int countCSVFileRows(const std::string& fileName)
{
    ifstream ifs0(fileName);
    if (!ifs0.is_open())
        throw runtime_error("unable to open input file for counting lines.\n");
    string str;
    int count = 0;
    while(getline(ifs0,str)){
        ++count;
    }
    return count;
}
