#pragma once
#include <potentials/potential.h>
#include <neighbourlists.h>

class LennardJones : public Potential {
private:
    double m_sigma;
    double m_epsilon;
public:
    LennardJones(double sigma, double epsilon);
    ~LennardJones() {}
    void calculateForces(System *system);
    void calculateForcesWithCells(System *system, Neighbourlists *neighbourlist);
};
