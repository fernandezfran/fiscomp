set encoding utf8

set xlabel "t"
set ylabel "Presión/N"
set grid
set key r b

p '../../../05-Molecular_dynamics/02/a-compare-bd/thermo.dat' u 1:($6/500) t 'Dinámica Molecular', \
  'thermo.dat' u 2:($4/512) t 'Dinámica Browniana'

set terminal postscript enhanced color eps 20
set output "pres.eps"
rep
set terminal qt

