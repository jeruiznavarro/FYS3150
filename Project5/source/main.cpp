#include <iostream>
#include <math/random.h>
#include <potentials/lennardjones.h>
#include <integrators/velocityverlet.h>
#include <integrators/eulercromer.h>
#include <system.h>
#include <statisticssampler.h>
#include <atom.h>
#include <io.h>
#include <unitconverter.h>
#include <thermostat.h>
#include <neighbourlists.h>

using namespace std;

int main() {
    double dt = UnitConverter::timeFromSI(1e-14);

    System system;
    VelocityVerlet *integrator = new VelocityVerlet();
    //EulerCromer *integrator = new EulerCromer();
    system.setIntegrator(integrator);
    system.createFCCLattice(5, UnitConverter::lengthFromAngstroms(5.26), 687.5);
    system.removeMomentum();
    LennardJones *potential = new LennardJones(3.405, 1.);
    system.setPotential(potential);
    StatisticsSampler statisticsSampler;
    Thermostat thermostat(&statisticsSampler, &system);
    thermostat.setThermostat(false);
    Neighbourlists neighbourlist(&system);

    IO *output = new IO();
    output->open((char *)"output.dat");
    IO *movie = new IO();
    movie->open((char *)"movie.xyz");
    movie->saveState(system);
    for(int timestep=0; timestep<1000; timestep++) {
        system.step(dt, thermostat, neighbourlist);
        statisticsSampler.sample(system);
        movie->saveState(system);
        if(timestep % 10 == 0) {
            output->saveOutput(system, statisticsSampler);
        }
        if(timestep % 100 == 0) {
            cout << timestep << endl;
        }
    }
    statisticsSampler.getEnergyVariance();
    cout << "The standard deviation of the energy is:" << statisticsSampler.energyTotalVariance << endl;
    movie->close();
    output->close();

    return 0;
}



//cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
//cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
//cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
//cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
//cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;
//cout << "One unit of pressure is " << UnitConverter::pressureToSI(1.0) << " Pa" << endl;
