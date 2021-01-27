#tiempo vs N
set encoding utf8

set xlabel "r"
set ylabel "g(r)"
set grid
set key r t

p 'gder08.dat' u 1:2 w lp pt 5 dt 2 t '{/Symbol r} = 0.8', \
  'gder12.dat' u 1:2 w lp pt 7 dt 2 t '{/Symbol r} = 1.2'

set terminal postscript enhanced color eps 20
set output "gder.eps"
rep
set terminal qt
