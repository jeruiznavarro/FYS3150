#include <system.h>
#include <math.h>
#include <math/vec3.h>
#include <integrators/integrator.h>
#include <potentials/potential.h>
#include <unitconverter.h>
#include <statisticssampler.h>
#include <thermostat.h>
#include <neighbourlists.h>

System::System() :
    m_potential(0),
    m_integrator(0),
    m_currentTime(0),
    m_steps(0)
{

}

System::~System() {
    delete m_potential;
    delete m_integrator;
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
    for(unsigned int i=0; i<m_atoms.size(); i++) {
        Atom *atom = m_atoms[i];
        for(unsigned int j=0; j<3; j++) {
            if(atom->position[j] < 0.) {
                atom->position[j] += systemSize()[j];
                atom->unwrappedPosition[j] -= systemSize()[j];
            }
            if(atom->position[j] >= systemSize()[j]) {
                atom->position[j] -= systemSize()[j];
                atom->unwrappedPosition[j] += systemSize()[j];
            }
        }
    }
}

void System::removeMomentum() {
    vec3 sum;
    for(unsigned int i=0; i<m_atoms.size(); i++) {
        Atom *atom = m_atoms[i];
        for (unsigned int j=0; j<3; j++) {
            sum[j] += atom->velocity[j];
        }
    }
    for (unsigned int j=0; j<3; j++) {
        sum[j] /= m_atoms.size();
    }
    for(unsigned int i=0; i<m_atoms.size(); i++) {
        Atom *atom = m_atoms[i];
        for (unsigned int j=0; j<3; j++) {
            atom->velocity[j] -= sum[j];
        }
    }
}

void System::resetForcesOnAllAtoms() {
    for(unsigned int i=0; i<m_atoms.size(); i++) {
        Atom *atom = m_atoms[i];
        atom->force.setToZero();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    for(int i=0; i<numberOfUnitCellsEachDimension; i++) {
        for(int j=0; j<numberOfUnitCellsEachDimension; j++) {
            for(int k=0; k<numberOfUnitCellsEachDimension; k++) {
                Atom *atom0 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom0->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(temperature));
                atom0->position[0] = i*latticeConstant;
                atom0->position[1] = j*latticeConstant;
                atom0->position[2] = k*latticeConstant;
                atom0->initialPosition[0] = atom0->position[0];
                atom0->initialPosition[1] = atom0->position[1];
                atom0->initialPosition[2] = atom0->position[2];
                atom0->unwrappedPosition[0] = atom0->position[0];
                atom0->unwrappedPosition[1] = atom0->position[1];
                atom0->unwrappedPosition[2] = atom0->position[2];
                atoms().push_back(atom0);
                Atom *atom1 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom1->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(temperature));
                atom1->position[0] = (i+0.5)*latticeConstant;
                atom1->position[1] = (j+0.5)*latticeConstant;
                atom1->position[2] = k*latticeConstant;
                atom1->initialPosition[0] = atom1->position[0];
                atom1->initialPosition[1] = atom1->position[1];
                atom1->initialPosition[2] = atom1->position[2];
                atom1->unwrappedPosition[0] = atom1->position[0];
                atom1->unwrappedPosition[1] = atom1->position[1];
                atom1->unwrappedPosition[2] = atom1->position[2];
                atoms().push_back(atom1);
                Atom *atom2 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom2->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(temperature));
                atom2->position[0] = i*latticeConstant;
                atom2->position[1] = (j+0.5)*latticeConstant;
                atom2->position[2] = (k+0.5)*latticeConstant;
                atom2->initialPosition[0] = atom2->position[0];
                atom2->initialPosition[1] = atom2->position[1];
                atom2->initialPosition[2] = atom2->position[2];
                atom2->unwrappedPosition[0] = atom2->position[0];
                atom2->unwrappedPosition[1] = atom2->position[1];
                atom2->unwrappedPosition[2] = atom2->position[2];
                atoms().push_back(atom2);
                Atom *atom3 = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                atom3->resetVelocityMaxwellian(UnitConverter::temperatureFromSI(temperature));
                atom3->position[0] = (i+0.5)*latticeConstant;
                atom3->position[1] = j*latticeConstant;
                atom3->position[2] = (k+0.5)*latticeConstant;
                atom3->initialPosition[0] = atom3->position[0];
                atom3->initialPosition[1] = atom3->position[1];
                atom3->initialPosition[2] = atom3->position[2];
                atom3->unwrappedPosition[0] = atom3->position[0];
                atom3->unwrappedPosition[1] = atom3->position[1];
                atom3->unwrappedPosition[2] = atom3->position[2];
                atoms().push_back(atom3);
            }
        }
    }
    m_dimensionSize = numberOfUnitCellsEachDimension*latticeConstant;
    setSystemSize(vec3(m_dimensionSize, m_dimensionSize, m_dimensionSize));
    m_systemVolume = m_dimensionSize*m_dimensionSize*m_dimensionSize;
}

void System::calculateForces(Neighbourlists neighbourlist) {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
    //m_potential->calculateForcesWithCells(this, &neighbourlist);
}

void System::step(double dt, Thermostat thermostat, Neighbourlists neighbourlist) {
    m_integrator->integrate(this, &thermostat, &neighbourlist, dt);
    m_steps++;
    m_currentTime += dt;
}
