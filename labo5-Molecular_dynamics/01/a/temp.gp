set encoding utf8

set xlabel 't'
set ylabel 'T'
set grid
set key r b

p 'thermo.dat' u 1:5 w lp pt 7 dt 2 t 'Temperatura instantanea', \


set terminal postscript enhanced color eps 20
set output "temp.eps"
rep
set terminal qt

