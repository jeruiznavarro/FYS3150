#pragma once
#include <integrators/integrator.h>
#include <neighbourlists.h>
#include <math/vec3.h>

class EulerCromer : public Integrator {
private:
    bool m_firstStep;
public:
    EulerCromer() {}
    ~EulerCromer() {}
    void integrate(System* system, Thermostat *thermostat, Neighbourlists *neighbourlist, double dt);
};
