#!bin/bash
#
#script para compilar el programa de integración numérica
#y obtener el gráfico de gnuplot.
#

gfortran -c precision.f90
gfortran -c -O3 gaussmod.f90
gfortran -c -O3 integrales.f90

gfortran -O3 main.f90 precision.o gaussmod.o integrales.o && ./a.out

gnuplot e_vs_n.gp
