set encoding utf8

#set title 'T = 2.0'
set xlabel 'MC step'
set ylabel 'E/N, M/N'
set grid
set key r c

set xrange [0:150]
p 'EM_2.00-mess.dat' u 1:3 w lp pt 7 dt 2 t 'Magnetización T_{0} desorden', \
  'EM_2.00-order.dat' u 1:3 w lp pt 9 dt 2 t 'Magnetización T_{0} orden', \
  'EM_2.00-mess.dat' u 1:2 w lp pt 11 dt 2 t 'Energía T_{0} desorden', \
  'EM_2.00-order.dat' u 1:2 w lp pt 13 dt 2 t 'Energía T_{0} orden'

set terminal postscript enhanced color eps 18
set output "EM_2.00.eps"
rep
set terminal qt
