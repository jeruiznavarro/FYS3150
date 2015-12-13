#include "thermostat.h"
#include "system.h"
#include "statisticssampler.h"
#include <cmath>

Thermostat::Thermostat(StatisticsSampler *sampler, System *system) {
    this->sampler = sampler;
    this->system = system;
}

void Thermostat::berendsen(double bathTemperature, double tau, double dt) {
    for (unsigned int i=0; i<system->atoms().size(); i++) {
        Atom *atom = system->atoms()[i];
        atom->velocity = atom->velocity*sqrt(1.+dt*(bathTemperature/sampler->temperature-1.)/tau);
    }
}

void Thermostat::setThermostat(bool value) {
    active = value;
}
void Thermostat::switchThermostat() {
    active = !active;
}
