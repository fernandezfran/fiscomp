set encoding utf8

set ylabel "{/Symbol r}"
set xlabel "t"
#set xrange [0:100]
#set yrange [1.05:2.15]
#set ytics 1.1,0.2,2.2
set grid ytics
set key box r t

p 'thermo.6.001000.dat' u 1:9 w p pt 5 t '{/Symbol bt}_P = 1000', \
  'thermo.6.010000.dat' u 1:9 w p pt 7 t '{/Symbol bt}_P = 10000', \
  'thermo.6.100000.dat' u 1:9 w p pt 9 t '{/Symbol bt}_P = 100000'

set term postscript enhanced color eps 22
set output "rhoterm.6.eps"
rep
set terminal qt
