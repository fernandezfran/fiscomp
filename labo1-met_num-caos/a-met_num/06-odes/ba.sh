#!/bin/bash
#
#script de bash con como compilé el programa y con los graficos de gnuplot
#el archivo in.metodo se va modificando para cambiar el algoritmo y poder
#medir cuanto tiempo se demora cada uno.

gfortran -c precision.f90
gfortran -c -O3 odes.f90

gfortran -O3 precision.o odes.o main_ab.f90 -o a-b.out

printf 1 > in.metodo
time ./a-b.out

printf 2 > in.metodo
time ./a-b.out

printf 3 > in.metodo
time ./a-b.out

gfortran -O3 precision.o odes.o main_c.f90 -o c.out && ./c.out

gnuplot energia.gp
python3 sort.py
gnuplot error_global.gp 
gnuplot error_local.gp
gnuplot x_v_t.gp
