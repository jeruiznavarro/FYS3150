#include <io.h>
#include <system.h>
#include <atom.h>
#include <unitconverter.h>
#include <statisticssampler.h>
#include <potentials/lennardjones.h>
#include <cstdlib>
#include <iomanip>
#include <iostream>


using namespace std;
using std::endl;
using std::cout;

IO::IO() {
    firstTime = true;
}

IO::~IO() {
    close();
}

void IO::open(char *filename) {
    if(file.is_open()) {
        std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
        exit(1);
    }

    file.open(filename);
}

void IO::close() {
    if(file.is_open()) {
        file.close();
    }
}

void IO::saveState(System &system) {
    file << system.atoms().size() << endl;
    file << "#This is an optional comment." << endl;
    for(unsigned int n=0; n<system.atoms().size(); n++) {
        Atom *atom = system.atoms()[n];
        file << "Ar " << UnitConverter::lengthToAngstroms(atom->position.x()) << " " << UnitConverter::lengthToAngstroms(atom->position.y()) << " " << UnitConverter::lengthToAngstroms(atom->position.z()) << endl;
    }
}
void IO::saveOutput(System &system, StatisticsSampler &statisticsSampler) {
    if(firstTime) {
        file << "# Time, total energy, potential energy, kinetic energy, temperature, pressure & MSD/6" << endl;
        firstTime = false;
    }
    file << std::setprecision(15) << UnitConverter::timeToSI(system.currentTime()) << " " << UnitConverter::energyToEv(statisticsSampler.potentialEnergy+statisticsSampler.kineticEnergy) << " " << UnitConverter::energyToEv(statisticsSampler.potentialEnergy) << " " << UnitConverter::energyToEv(statisticsSampler.kineticEnergy) << " " << UnitConverter::temperatureToSI(statisticsSampler.temperature) << " " << UnitConverter::pressureToSI(statisticsSampler.pressure) << " " << UnitConverter::areaToSI(statisticsSampler.meanSquareDisplacement)/(6*system.atoms().size()) << endl;
}

//void IO::saveToBinary(string *filename, System *system) {
//    ofstream file(filename, ios::in | ios::binary);
//    if (!file.is_open()) {
//        cerr << "Could not open file " << filename << ". Aborting!" << endl;
//        exit (1);
//    }
//    int numberOfPhseSpaceCoordinates = 6*system->atoms().size();
//    double *phaseSpace  = new double [numberOfPhseSpaceCoordinates];
//    vec3 systemSize = system->systemSize();
//    int phaseSpaceCounter = 0;
//    for(unsigned int i=0; i<system->atoms().size(); i++) {
//        Atom *atom = system->atoms()[i];
//        phaseSpace[phaseSpaceCounter++] = atom->position[0];
//        phaseSpace[phaseSpaceCounter++] = atom->position[1];
//        phaseSpace[phaseSpaceCounter++] = atom->position[2];
//        phaseSpace[phaseSpaceCounter++] = atom->velocity[0];
//        phaseSpace[phaseSpaceCounter++] = atom->velocity[1];
//        phaseSpace[phaseSpaceCounter++] = atom->velocity[2];
//    }
//    int numberofAtoms = system->atoms().size();
//    file.write(reinterpret_cast<char*>(numberofAtoms), sizeof(int));
//    file.write(reinterpret_cast<char*>(phaseSpace), numberOfPhseSpaceCoordinates*sizeof(double));
//    file.write(reinterpret_cast<char*>(&systemSize[0]), 3*sizeof(double));
//    file.close();
//    delete phaseSpace;
//}

//void IO::loadFromBinary(string *filename, System *system) {
//    ofstream file(filename, ios::in | ios::binary);
//    if (!file.is_open()) {
//        cerr << "Could not open file " << filename << " . Aborting!" << endl;
//        exit (1);
//    }
//    int numberOfAtoms ;
//    file.read(reinterpret_cast<char*>(&numberOfAtoms), sizeof(int));
//    double *phaseSpace  = new double [6*numberOfAtoms];
//    file.read(reinterpret_cast<char*>(phaseSpace), 6*numberOfAtoms*sizeof(double));
//    vec3 systemSize;
//    file.read(reinterpret_cast<char*>(&systemSize[0], 3*sizeof(double)));
//    system->setSystemSize(systemSize);
//    int phaseSpaceCounter = 0 ;
//    double mass = 39.948;
//    system->atoms().clear();
//    for (int i =0; i <numberOfAtoms; i++) {
//        Atom *atom = new Atom(mass);
//        atom->position = vec3(phaseSpace[phaseSpaceCounter++], phaseSpace[phaseSpaceCounter++], phaseSpace[phaseSpaceCounter++]);
//        atom->velocity = vec3(phaseSpace[phaseSpaceCounter++], phaseSpace[phaseSpaceCounter++], phaseSpace[phaseSpaceCounter++]);
//        system->atoms().push_back(atom);
//    }
//    file.close();
//    delete phaseSpace;
//}
