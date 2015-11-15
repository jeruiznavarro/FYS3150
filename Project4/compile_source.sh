#!/bin/bash
gfortran  -Wall -Wextra -Wtabs -mcmodel=large -fPIC -g -fcheck=all -fbacktrace -ffree-line-length-0 ising_2dim.f90 f90lib.f90 -o 4_debug.exe #debug
gfortran  -Wall -Wextra -Wtabs -mcmodel=large -fPIC -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-0 ising_2dim.f90 f90lib.f90 -o 4.exe #production run
#rm -f *.mod *.o