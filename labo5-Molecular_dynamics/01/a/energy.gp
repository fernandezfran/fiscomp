set encoding utf8

set xlabel 't'
set ylabel 'E'
set grid
set key at 10,-100 


p 'thermo.dat' u 1:2 w lp pt 9 dt 2 t 'Energía cinética', \
  ''           u 1:3 w lp pt 7 dt 2 t 'Energía potencial', \
  ''           u 1:4 w lp pt 5 dt 2 t 'Energía total'


set terminal postscript enhanced color eps 20
set output "energy.eps"
rep
set terminal qt
