#pragma once
#include <vector>
#include <atom.h>
#include <cmath>
using std::vector;

class System;
class Neighbourlists {
private:
    vector<vector<vector<vector<Atom*> > > > m_cells;
    double m_cutOffSize;
    double m_cellSize;
    unsigned int m_cellsPerDimension;
public:
    Neighbourlists(System *system);
    void listCells(System *system);
    void clearCells();
    void updateCells(System *system);
    vector<Atom*> &cellWithIndex(unsigned int i, unsigned int j, unsigned int k) { return m_cells[i][j][k]; }
    int getCellOccupancy(unsigned int i, unsigned int j, unsigned int k) { return m_cells[i][j][k].size(); }
    double getCellSize() { return m_cellSize; }
    double getCutOffSize() { return m_cutOffSize; }
    double getCellsPerDimension() { return m_cellsPerDimension; }
};
