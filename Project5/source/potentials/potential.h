#pragma once
#include <string>
#include <vector>
#include <system.h>
#include <neighbourlists.h>

class Potential {
protected:
    double m_potentialEnergy;
    double m_pressure;
public:
    Potential();
    virtual ~Potential() {}
    virtual void calculateForces(System *system) = 0;
    virtual void calculateForcesWithCells(System *system, Neighbourlists *neighbourlist) = 0;
    double potentialEnergy() { return m_potentialEnergy; }
    double pressureFromPotential() { return m_pressure; }
};
