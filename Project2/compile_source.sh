#!/bin/bash
gfortran f90lib.f90 project2_1.f95 -o 1.exe
gfortran project2_2.f95 -o 2.exe
gfortran project2_test.f95 -o test.exe
rm -f *.mod *.o