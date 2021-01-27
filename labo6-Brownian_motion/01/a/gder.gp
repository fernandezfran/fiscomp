#tiempo vs N
set encoding utf8

set xlabel "r"
set ylabel "g(r)"
set grid
set key r t

p '../../../05-Molecular_dynamics/02/a-compare-bd/gder.dat' u 1:2 w lp pt 5 ps 1.2 dt 2 t 'Dinámica Molecular', \
  'gder.dat' u 1:2 w lp pt 7 ps 0.8 dt 2 t 'Dinámica Browniana'

set terminal postscript enhanced color eps 20
set output "gder.eps"
rep
set terminal qt

