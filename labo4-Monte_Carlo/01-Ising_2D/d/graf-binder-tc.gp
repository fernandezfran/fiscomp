set encoding utf8

set xlabel 'T/T_C'
set ylabel 'Cumulante de Binder'
set grid
set key r t 
set xrange [0.97:1.03]
set yrange [0.55:0.67]

set arrow from 0.998,graph(0,0) to 0.998,graph(1,1) nohead lc 6 

p 'tul10.dat' u 1:(10*10*$2) w lp pt 11 dt 2 t 'L=10', \
  'tul20.dat' u 1:(20*20*$2) w lp pt  9 dt 2 t 'L=20', \
  'tul40.dat' u 1:(40*40*$2) w lp pt  7 dt 2 t 'L=40', \
  'tul128.dat' u 1:(128*128*$2) w lp pt 5 dt 2 t 'L=128', \
  0 w l lc 6 t 'T = 0.998*T_C'


set terminal postscript enhanced color eps 18
set output "graf-binder-tc.eps"
rep
set terminal qt
