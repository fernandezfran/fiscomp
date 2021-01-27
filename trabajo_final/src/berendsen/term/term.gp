set encoding utf8

set ylabel "T"
set xlabel "t"
set xrange [0:100]
set yrange [1.05:2.15]
set ytics 1.1,0.2,2.2
set grid ytics
set key box r b

p 'thermo0010.dat' u 1:5 w p pt 5 t '{/Symbol t}_T =10', \
  'thermo0100.dat' u 1:5 w p pt 7 t '{/Symbol t}_T = 100', \
  'thermo1000.dat' u 1:5 w p pt 9 t '{/Symbol t}_T = 1000'

set term postscript enhanced color eps 22
set output "term.eps"
rep
set terminal qt
