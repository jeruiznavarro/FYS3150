#!/bin/bash
gfortran  -Wall -Wextra -Wtabs -mcmodel=large -fPIC -g -fcheck=all -fbacktrace -ffree-line-length-0 f90lib.f90 gauss-laguerre.f90 dranxor2.f95 project3.f08 -o 3_debug.exe #debug
gfortran  -Wall -Wextra -Wtabs -mcmodel=large -fPIC -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-0 f90lib.f90 gauss-laguerre.f90 dranxor2.f95 project3.f08 -o 3.exe #production run
gfortran  -Wall -Wextra -Wtabs -mcmodel=large -fPIC -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-0 f90lib.f90 test_legendre.f08 -o test_legendre.exe #test
gfortran  -Wall -Wextra -Wtabs -mcmodel=large -fPIC -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-0 gauss-laguerre.f90 test_laguerre.f08 -o test_laguerre.exe #test
gfortran  -Wall -Wextra -Wtabs -mcmodel=large -fPIC -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-0 dranxor2.f95 test_mc.f08 -o test_mc.exe #test
rm -f *.mod *.o