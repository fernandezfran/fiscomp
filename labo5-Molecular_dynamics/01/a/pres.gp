set encoding utf8

set xlabel 't'
set ylabel 'P'
set grid
set key r b

p 'thermo.dat' u 1:6 w lp pt 7 dt 2 t 'Presión instantanea', \


set terminal postscript enhanced color eps 20
set output "pres.eps"
rep
set terminal qt
