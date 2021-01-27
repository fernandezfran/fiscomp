#Factor de estructura para dos densidades
set encoding utf8

set xlabel "t"
set ylabel "S(k)"
set grid
set format x "%g"
set format y "10^{%T}"
set log y
set key l b 

p 'thermo08.dat' u 1:2 w lp pt 5 dt 2 t '{/Symbol r} = 0.8', \
  'thermo12.dat' u 1:2 w lp pt 7 dt 2 t '{/Symbol r} = 1.2'

set terminal postscript enhanced color eps 20
set output "sk.eps"
rep
set terminal qt

