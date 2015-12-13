#include <statisticssampler.h>
#include <system.h>
#include <potentials/potential.h>
#include <math.h>

StatisticsSampler::StatisticsSampler() {
    energySumCounter = 0;
}

StatisticsSampler::~StatisticsSampler() {

}

void StatisticsSampler::sample(System &system) {
    sampleKineticEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    samplePressure(system);
    sampleHeatCapacity(system);
    sampleDiffusionConstant(system);
    potentialEnergy = system.potential()->potentialEnergy();
    sampleEnergyVariance();
}

void StatisticsSampler::sampleKineticEnergy(System &system) {
    kineticEnergy = 0;
    for (unsigned int i=0; i<system.atoms().size(); i++) {
        Atom *atom = system.atoms()[i];
        kineticEnergy += 0.5*atom->mass()*atom->velocity.lengthSquared();
    }
    energySum += kineticEnergy;
    energySumSquared += kineticEnergy*kineticEnergy;
    energySumCounter += 1;
}

void StatisticsSampler::sampleTemperature(System &system) {
    temperature = (2.*kineticEnergy)/(3.*system.atoms().size());
}

void StatisticsSampler::sampleDensity(System &system) {
    density = system.atoms().size()/system.systemVolume();
}

void StatisticsSampler::samplePressure(System &system) {
    pressure = system.potential()->pressureFromPotential()+density*temperature;
}

void StatisticsSampler::sampleHeatCapacity(System &system) {
    heatCapacity = 3./(2.-4.*system.atoms().size()*((energySumSquared-energySum*energySum/energySumCounter)/energySumCounter)/(3.*temperature*temperature));
}

void StatisticsSampler::sampleDiffusionConstant(System &system) {
    for (unsigned int i=0; i<system.atoms().size(); i++) {
        Atom *atom = system.atoms()[i];
        meanSquareDisplacement += (atom->unwrappedPosition-atom->initialPosition).lengthSquared();
    }
}

void StatisticsSampler::sampleEnergyVariance() {
    energyTotalSum += kineticEnergy+potentialEnergy;
    energyTotalSumSquared += (kineticEnergy+potentialEnergy)*(kineticEnergy+potentialEnergy);
}

void StatisticsSampler::getEnergyVariance() {
    energyTotalVariance = sqrt(energyTotalSumSquared/energySumCounter-(energyTotalSum*energyTotalSum)/(energySumCounter*energySumCounter));
}
