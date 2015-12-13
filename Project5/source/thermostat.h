#pragma once
#include <statisticssampler.h>
#include <system.h>

class Thermostat {
private:
    StatisticsSampler *sampler;
    System* system;

public:
    Thermostat(StatisticsSampler* sampler, System* system);
    void berendsen(double bathTemperature, double tau, double dt);
    void setThermostat(bool value);
    void switchThermostat();
    bool active;
};

