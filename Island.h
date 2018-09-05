#ifndef ISLAND_H
#define ISLAND_H

class Island {        // class for ONE island within archipelago
public:
    explicit Island(const int k) : mIK{k} {assert(k >= 0);} // island constructor based on island-wide K

    int get_carrying_capacity() { return mIK; }
    // int sizeIsl() const {return static_cast<int>(mvIsland.size());}  // returns size of island vector
    int specAlive() const {return static_cast<int>(mvIslSpecAlive.size());} // returns number of species alive
    double returnLogGrowth() { return 1 - static_cast<double>(mvIslSpecAlive.size()) / mIK;}   // returns the logistic growth term (1-n/K)

    int findPos(const int &ID) const;    // find the position of certain species (input) in IslandPhylo vector

    const int createNewID();       // returns new species ID and maxID += 1

    double calculateIslRates(const vector<double> &, const int &, const int &, const double &);
    // initialise/calculate rates and store them in EventRates vector
    // gam_i, gam_m, lamb_cl, lamb_al, mu_l
    // per island -> doesn't include global rates !!!
    const double extractSumIslRate() const noexcept;      // return the per-island rates vector

    vector<int> sampleLocalEvent(mt19937_64, const int&);   // in case a local event is drawn, sample island, event and species
    // it happens to

    void immigrate(const int&, const double& BirthT);                   // mainland species immigrates to that island
    int migrate(const int &, vector<double> &, const double &, mt19937_64);                     // island species migrates to other island
    void speciateClado(const int &);               // island species cladogenetically speciates
    void speciateAna(const int &);                 // island species anagenetically speciates
    void goExtinct(const int &);                   // island species goes extinct

    Species returnSpecies(const int &iPos) { return mvIsland[iPos]; }   // returns specific species from species vector
    void pushbackSp(const Species &spNew) { mvIsland.push_back(spNew); }    // adds new species to species vector

    void addIsland(const Island &);     // add island to THIS islandd

    void printIsland();                 // prints island vector of species to the screen

    const vector<Species>& returnIsland() const { return mvIsland; }    // return island vector of species
    const vector<int>& returnIslSpecAlive() const { return mvIslSpecAlive; }  // return extant species vector
private:
    vector<Species> mvIsland;           // phylogeny vector of species on island
    vector<int> mvIslSpecAlive;            // vector of species' IDs of alive species (.size() = n species on island)
    vector<double> mvLocalRates;        // vector of rates for events PER ISLAND (5 rates: gamI, gamM, lambC, lambA, mu)
    int mIK; // Carrying capacity (should be const one day)
    // for now: mIk = mAK / iNumIslands
};

#endif // ISLAND_H
