#pragma once
#include <integrators/integrator.h>
#include <neighbourlists.h>

class VelocityVerlet : public Integrator {
private:
    void firstKick(System *system, double dt);
    void halfKick(System *system, double dt);
    void move(System *system, Neighbourlists *neighbourlist, double dt);
    bool m_firstStep;
public:
    VelocityVerlet();
    ~VelocityVerlet();
    void integrate(System *system, Thermostat *thermostat, Neighbourlists *neighbourlist, double dt);
};
