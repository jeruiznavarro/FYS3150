#pragma once
#include <fstream>
class System;
class StatisticsSampler;
using std::ofstream;
using namespace std;

class IO {
private:
    ofstream file;
    bool firstTime;
public:
    IO();
    ~IO();
    void saveState(System &system);
    void saveOutput(System &system, StatisticsSampler &statisticsSampler);
    void open(char *filename);
    void close();
    //void saveToBinary(string *filename, System *system);
    //void loadFromBinary(string *filename, System *system);
};
