#!bin/bash

gfortran -c precision.f90
gfortran -c -O3 dft_mod.f90

gfortran -O3 main.f90 precision.o dft_mod.o -lfftw3

./a.out > fftw_vs_dft.dat

gnuplot t_vs_N.gp
