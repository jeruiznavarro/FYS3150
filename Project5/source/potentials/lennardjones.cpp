#include <iostream>
#include <potentials/lennardjones.h>
#include <cmath>
#include <math/vec3.h>
#include <io.h>
#include <system.h>
#include <unitconverter.h>
#include <neighbourlists.h>
#include <vector>
#include <algorithm>

using namespace std;

LennardJones::LennardJones(double sigma, double epsilon) {
    m_sigma = sigma;
    m_epsilon = epsilon;
}

void LennardJones::calculateForces(System *system) {
    int counter = 0;
    m_potentialEnergy = 0;
    for(unsigned int i=0; i<system->atoms().size(); i++) {
        Atom *atom0 = system->atoms()[i];
        for(unsigned int j=i+1; j<system->atoms().size(); j++) {
            Atom *atom1 = system->atoms()[j];
            vec3 deltaRVector = atom0->position-atom1->position;
            for(unsigned int k=0; k<3; k++) {
                if(deltaRVector[k] > system->systemSize()[k]*0.5) deltaRVector[k] -= system->systemSize()[k];
                else if(deltaRVector[k] < -system->systemSize()[k]*0.5) deltaRVector[k] += system->systemSize()[k];
            }
            double distance = deltaRVector.length();
            double force = 24.*m_epsilon*pow(m_sigma/distance,6)*(2*pow(m_sigma/distance,6)-1.)/(distance*distance);
            m_pressure -= force*distance;
            m_potentialEnergy += 4.*m_epsilon*pow(m_sigma/distance,6)*(pow(m_sigma/distance,6)-1.);
            counter++;
            atom0->force.addAndMultiply(deltaRVector, force);
            atom1->force.addAndMultiply(deltaRVector, -force);
        }
    }
    m_pressure /= (3.*counter*system->systemVolume());
    counter = 0;
}

void LennardJones::calculateForcesWithCells(System *system, Neighbourlists *neighbourlist) {
    vec3 displacementVector;
    vec3 distanceVector;
    int ll, mm, nn;
    int cells = neighbourlist->getCellsPerDimension(), counter = 0;
    double size = system->dimensionSize();
    m_potentialEnergy = 0.;
    m_pressure = 0.;
    for (int i=0; i<cells; i++) {
        for (int j=0; j<cells; j++) {
            for (int k=0; k<cells; k++) {
                vector<Atom*> cell0 = neighbourlist->cellWithIndex(i, j, k);
                for (int ii=0; ii<neighbourlist->getCellOccupancy(i, j, k); ii++) {
                    Atom *atom0 = cell0[ii];
                    for (int l=i-1; l<i+2; l++) {
                        if ((i == 0) && (l == i-1)) {
                            ll = cells-1;
                            displacementVector[0] = -size;
                        }
                        else if ((i == cells-1) && (l == i+1)) {
                            ll = 0;
                            displacementVector[0] = size;
                        }
                        else {
                            ll = l;
                            displacementVector[0] = 0.;
                        }
                        for (int m=j-1; m<j+2; m++) {
                            if ((j == 0) && (m == j-1)) {
                                mm = cells-1;
                                displacementVector[1] = -size;
                            }
                            else if ((j == cells-1) && (m == j+1)) {
                                mm = 0;
                                displacementVector[1] = size;
                            }
                            else {
                                mm = m;
                                displacementVector[1] = 0.;
                            }
                            for (int n=k-1; n<k+2; n++) {
                                if ((k == 0) && (n == k-1)) {
                                    nn = cells-1;
                                    displacementVector[2] = -size;
                                }
                                else if ((k == cells-1) && (n == k+1)) {
                                    nn = 0;
                                    displacementVector[2] = size;
                                }
                                else {
                                    nn = n;
                                    displacementVector[2] = 0.;
                                }
                                vector<Atom*> cell1 = neighbourlist->cellWithIndex(ll, mm, nn);
                                for (int jj=0; jj<neighbourlist->getCellOccupancy(ll, mm, nn); jj++) {
                                    if ((ii == jj) && (i == ll) && (j == mm) && (k == nn)) continue;
                                    Atom *atom1 = cell1[jj];
                                    for (unsigned int kk=0; kk<3; kk++) {
                                        if (abs(atom0->position[kk]-atom1->position[kk]) < size*0.5) {
                                            displacementVector[kk] = 0.;
                                        }
                                        distanceVector[kk] = atom0->position[kk]-atom1->position[kk]-displacementVector[kk];
                                    }
                                    double distance = distanceVector.length();
                                    if (distance > neighbourlist->getCutOffSize()) continue;
                                    counter++;
                                    double force = 24.*m_epsilon*pow(m_sigma/distance,6)*(2*pow(m_sigma/distance,6)-1.)/(distance*distance);
                                    m_pressure -= force*distance;
                                    m_potentialEnergy += 4.*m_epsilon*pow(m_sigma/distance,6)*(pow(m_sigma/distance,6)-1.)+0.02;
                                    atom0->force.addAndMultiply(distanceVector, force);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    m_potentialEnergy /= 2.;
    m_pressure /= (6.*counter*system->systemVolume());
    counter = 0;
}
