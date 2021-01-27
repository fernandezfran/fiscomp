set encoding utf8

set xlabel 't'
set ylabel 'Fluctuaciones de la temperatura'
set grid


p 'thermo.dat' u 1:(abs($5-0.109559E+01)) w lp pt 7 dt 2 t '', \


set terminal postscript enhanced color eps 20
set output "temp.eps"
rep
set terminal qt

