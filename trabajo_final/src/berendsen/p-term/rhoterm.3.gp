set encoding utf8

set ylabel "{/Symbol r}"
set xlabel "t"
#set xrange [0:100]
#set yrange [1.05:2.15]
#set ytics 1.1,0.2,2.2
set grid ytics
set key box r b

p 'thermo.3.001000.dat' u 1:9 w p pt 5 t '{/Symbol bt}_P = 1000', \
  'thermo.3.010000.dat' u 1:9 w p pt 7 t '{/Symbol bt}_P = 10000', \
  'thermo.3.100000.dat' u 1:9 w p pt 9 t '{/Symbol bt}_P = 100000', \

set term postscript enhanced color eps 22
set output "rhoterm.3.eps"
rep
set terminal qt
