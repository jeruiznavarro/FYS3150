#include <integrators/velocityverlet.h>
#include <system.h>
#include <thermostat.h>
#include <neighbourlists.h>

VelocityVerlet::VelocityVerlet() {
    m_firstStep = true;
}

VelocityVerlet::~VelocityVerlet() {

}

void VelocityVerlet::halfKick(System *system, double dt) {
    for(unsigned int i=0; i<system->atoms().size(); i++) {
        Atom *atom = system->atoms()[i];
        atom->velocity.addAndMultiply(atom->force, dt*0.5/atom->mass());
    }
}

void VelocityVerlet::move(System *system, Neighbourlists *neighbourlist, double dt) {
    for(unsigned int i=0; i<system->atoms().size(); i++) {
        Atom *atom = system->atoms()[i];
        atom->position.addAndMultiply(atom->velocity, dt);
    }
    system->applyPeriodicBoundaryConditions();
    if(m_firstStep) {
        neighbourlist->updateCells(system);
        m_firstStep = false;
    }
    else {
        neighbourlist->listCells(system);
    }
}

void VelocityVerlet::integrate(System *system, Thermostat *thermostat, Neighbourlists *neighbourlist, double dt) {
    if(m_firstStep) {
        neighbourlist->listCells(system);
        system->calculateForces(*neighbourlist);
    }
    halfKick(system, dt);
    move(system, neighbourlist, dt);
    system->calculateForces(*neighbourlist);
    halfKick(system, dt);
    if(thermostat->active) {
        thermostat->berendsen(10., dt, dt);
    }
}
