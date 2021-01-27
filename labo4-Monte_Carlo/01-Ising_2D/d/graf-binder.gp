set encoding utf8

set xlabel 'T/T_C'
set ylabel 'Cumulante de Binder'
set grid
set key r t 

p 'tul10.dat' u 1:(10*10*$2) w lp pt 11 dt 2 t 'L=10', \
  'tul20.dat' u 1:(20*20*$2) w lp pt  9 dt 2 t 'L=20', \
  'tul40.dat' u 1:(40*40*$2) w lp pt  7 dt 2 t 'L=40', \
  'tul128.dat' u 1:(128*128*$2) w lp pt 5 dt 2 t 'L=128'


set terminal postscript enhanced color eps 18
set output "graf-binder.eps"
rep
set terminal qt
