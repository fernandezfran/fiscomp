gfortran -c precision.f90

gfortran -O3 precision.o main.f90 -lfftw3

./a.out
