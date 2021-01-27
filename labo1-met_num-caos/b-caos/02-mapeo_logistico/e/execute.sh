gfortran -c precision.f90
gfortran -O3 main.f90 precision.o

./a.out
