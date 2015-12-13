#pragma once

class System;
class Thermostat;
class Neighbourlists;
class Integrator {
public:
    Integrator();
    virtual ~Integrator() { }
    virtual void integrate(System *system, Thermostat *thermostat, Neighbourlists *neighbourlist, double dt) = 0;
};
