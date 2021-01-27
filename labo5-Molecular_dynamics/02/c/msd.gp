#msd vs t
set encoding utf8

set xlabel "t"
set ylabel "MSD"
set grid
set format x "%g"
set format y "10^{%T}"
set log xy
set key l t

p 'msdvst-08-nve.dat' u 1:2 w lp pt 5 dt 2 t '{/Symbol r} = 0.8', \
  'msdvst-12-nve.dat' u 1:2 w lp pt 7 dt 2 t '{/Symbol r} = 1.2'

set terminal postscript enhanced color eps 20
set output "msd.eps"
rep
set terminal qt
