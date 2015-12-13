#pragma once
#include <system.h>

class StatisticsSampler {
public:
    double kineticEnergy;
    double potentialEnergy;
    double pressure;
    double density;
    double temperature;
    double heatCapacity;
    double energySum;
    double energySumSquared;
    double energyTotalSum;
    double energyTotalSumSquared;
    double energyTotalVariance;
    double meanSquareDisplacement;
    unsigned int energySumCounter;
    StatisticsSampler();
    ~StatisticsSampler();
    void sample(System &system);
    void sampleKineticEnergy(System &system);
    void sampleTemperature(System &system);
    void sampleDensity(System &system);
    void samplePressure(System &system);
    void sampleHeatCapacity(System &system);
    void sampleDiffusionConstant(System &system);
    void sampleEnergyVariance();
    void getEnergyVariance();
};
