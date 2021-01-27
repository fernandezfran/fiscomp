#!bin/bash
gfortran -c precision.f90
gfortran -c -O3 randomnum.f90
gfortran -c -O3 initial_cond.f90
gfortran -c -O3 calculate.f90
gfortran *.o -O3 -g -ftrapv main.f90
