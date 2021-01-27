gfortran -c precision.f90
gfortran -c odes.f90

gfortran main.f90 precision.o odes.o

printf 1 > tar
time ./a.out

printf 2 > tar
time ./a.out

printf 3 > tar
time ./a.out
