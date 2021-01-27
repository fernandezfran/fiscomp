set encoding utf8

set xlabel "t"
set ylabel "Energía potencial/N"
set grid
set key r t

p '../../../05-Molecular_dynamics/02/a-compare-bd/thermo.dat' u 1:($3/500) t 'Dinámica Molecular', \
  'thermo.dat' u 2:($3/512) t 'Dinámica Browniana'

set terminal postscript enhanced color eps 20
set output "epot.eps"
rep
set terminal qt

