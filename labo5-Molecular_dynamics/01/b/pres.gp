set encoding utf8

set xlabel 't'
set ylabel 'Fluctuaciones de la presión'
set grid


p 'thermo.dat' u 1:(abs($6-0.155018E+01)) w lp pt 7 dt 2 t '', \


set terminal postscript enhanced color eps 20
set output "pres.eps"
rep
set terminal qt
