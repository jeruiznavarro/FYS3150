#include <neighbourlists.h>
#include <system.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

Neighbourlists::Neighbourlists(System *system) {
    m_cutOffSize = 8.2273;
    m_cellsPerDimension = floor(system->dimensionSize()/m_cutOffSize);
    m_cellSize = system->dimensionSize()/m_cellsPerDimension;
    m_cells.resize(m_cellsPerDimension);
    for(unsigned int i=0; i<m_cellsPerDimension; i++) {
        m_cells[i].resize(m_cellsPerDimension);
        for(unsigned int j=0; j<m_cellsPerDimension; j++) {
            m_cells[i][j].resize(m_cellsPerDimension);
        }
    }
}

void Neighbourlists::listCells(System *system) {
    unsigned int index[3];
    for (unsigned int i=0; i<system->atoms().size(); i++) {
        Atom *atom = system->atoms()[i];
        for (int j=0; j<3; j++) {
             index[j] = floor(abs(atom->position[j])/m_cellSize);
        }
        m_cells[index[0]][index[1]][index[2]].push_back(atom);
    }
}

void Neighbourlists::clearCells() {
    for (unsigned int i=0; i<m_cellsPerDimension; i++) {
        for (unsigned int j=0; j<m_cellsPerDimension; j++) {
            for (unsigned int k=0; k<m_cellsPerDimension; k++) {
                m_cells[i][j][k].clear();
            }
        }
    }
}

void Neighbourlists::updateCells(System *system) {
    unsigned int a, c, newIndex[3], oldIndex[3];
    bool change;
    for (unsigned int i=0; i<m_cellsPerDimension; i++) {
        for (unsigned int j=0; j<m_cellsPerDimension; j++) {
            for (unsigned int k=0; k<m_cellsPerDimension; k++) {
                if (m_cells[i][j][k].size()<1) {
                    c++;
                }
            }
        }
    }
    for (unsigned int i=0; i<m_cellsPerDimension; i++) {
        oldIndex[0] = i;
        for (unsigned int j=0; j<m_cellsPerDimension; j++) {
            oldIndex[1] = j;
            for (unsigned int k=0; k<m_cellsPerDimension; k++) {
                oldIndex[2] = k;
                if (c == m_cellsPerDimension*m_cellsPerDimension*m_cellsPerDimension) {
                    listCells(system);
                }
                else {
                    for (unsigned int l=0; l<m_cells[i][j][k].size(); l++) {
                        Atom *atom = m_cells[i][j][k][l];
                        for (int m=0; m<3; m++) {
                            newIndex[m] = floor(abs(atom->position[m])/m_cellSize);
                            if (newIndex[m] != oldIndex[m]) {
                                change = true;
                            }
                        }
                        if (change) {
                            m_cells[newIndex[0]][newIndex[1]][newIndex[2]].push_back(atom);
                            m_cells[i][j][k].erase(m_cells[i][j][k].begin()+l);
                            change = false;
                        }
                    }
                }
                c = 0;
            }
        }
    }
    a  = system->atoms().size();
    a++;
}
