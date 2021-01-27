set encoding utf8

set ylabel "P(v)"
set xlabel "v"
set grid ytics
set key box
set yrange [0:.40]
set xrange [-5:5]

p (sqrt(1./2./2./pi))*exp(-x**2/(2.*2.)) lw 1.5 lc 8 t 'Maxwell-Boltzmann',\
  'histo0010.dat' u 1:2 pt 5 ps 1.2 t '{/Symbol t}_T = 10', \
  'histo0100.dat' u 1:2 pt 7 ps 1.0 t '{/Symbol t}_T = 100', \
  'histo1000.dat' u 1:2 pt 9 ps 0.8 t '{/Symbol t}_T = 1000'

set term postscript enhanced color eps 22
set output "veldist.eps"
rep
set terminal qt
