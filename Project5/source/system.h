#pragma once
#include <vector>
#include <atom.h>
#include <neighbourlists.h>
#include <math/vec3.h>

class Potential;
class Integrator;
class Thermostat;
using std::vector;
using CompPhys::vec3;

class System {
private:
    vec3 m_systemSize;
    vector<Atom*> m_atoms;
    Potential *m_potential;
    Integrator *m_integrator;
    double m_currentTime;
    double m_systemVolume;
    double m_dimensionSize;
    int m_steps;
public:
    System();
    ~System();
    void resetForcesOnAllAtoms();
    void createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature);
    void applyPeriodicBoundaryConditions();
    void removeMomentum();
    void calculateForces(Neighbourlists neighbourlist);
    void step(double dt, Thermostat thermostat, Neighbourlists neighbourlist);
    vector<Atom *> &atoms() { return m_atoms; }
    vec3 systemSize() { return m_systemSize; }
    void setSystemSize(vec3 systemSize) { m_systemSize = systemSize; }
    double systemVolume() { return m_systemVolume; }
    double dimensionSize() { return m_dimensionSize; }
    Potential *potential() { return m_potential; }
    void setPotential(Potential *potential) { m_potential = potential; }
    double currentTime() { return m_currentTime; }
    void setCurrentTime(double currentTime) { m_currentTime = currentTime; }
    Integrator *integrator() { return m_integrator; }
    void setIntegrator(Integrator *integrator) { m_integrator = integrator; }
    int steps() { return m_steps; }
    void setSteps(int steps) { m_steps = steps; }
};
