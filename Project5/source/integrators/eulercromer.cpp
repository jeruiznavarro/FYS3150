#include <integrators/eulercromer.h>
#include <system.h>
#include <thermostat.h>
#include <neighbourlists.h>
#include <math/vec3.h>

void EulerCromer::integrate(System* system, Thermostat *thermostat, Neighbourlists *neighbourlist, double dt) {
    if(m_firstStep) {
        neighbourlist->listCells(system);
        system->calculateForces(*neighbourlist);
    }
    for(unsigned int i=0; i<system->atoms().size(); i++) {
        Atom *atom = system->atoms()[i];
        //atom->velocity += atom->force*dt/atom->mass();
        //atom->position += atom->velocity*dt;
        //atom->velocity[0] += atom->force[0]*dt/atom->mass();
        //atom->velocity[1] += atom->force[1]*dt/atom->mass();
        //atom->velocity[2] += atom->force[2]*dt/atom->mass();
        //atom->position[0] += atom->velocity[0]*dt;
        //atom->position[1] += atom->velocity[1]*dt;
        //atom->position[2] += atom->velocity[2]*dt;
        atom->velocity.addAndMultiply(atom->force, dt/atom->mass());
        atom->position.addAndMultiply(atom->velocity, dt);
        atom->unwrappedPosition.addAndMultiply(atom->velocity, dt);
    }
    system->applyPeriodicBoundaryConditions();
    if(m_firstStep) {
        neighbourlist->updateCells(system);
        m_firstStep = false;
    }
    else {
        neighbourlist->listCells(system);
    }
    system->calculateForces(*neighbourlist);
    if(thermostat->active) {
        thermostat->berendsen(10., dt, dt);
    }
}
