set encoding utf8

set ylabel "T"
set xlabel "t"
set grid ytics
set key box r b

p 'thermo20.dat' u 1:5 w p pt 5 t '{/Symbol n} = 20', \
  'thermo10.dat' u 1:5 w p pt 7 t '{/Symbol n} = 10', \
  'thermo05.dat' u 1:5 w p pt 9 t '{/Symbol n} = 5', \
  'thermo01.dat' u 1:5 w p pt 11 t '{/Symbol n} = 1', \
  'thermo.1.dat' u 1:5 w p pt 13 t '{/Symbol n} = 0.1', \
  'thermo.01.dat' u 1:5 w p pt 15 t '{/Symbol n} = 0.01'

set term postscript enhanced color eps 22
set output "term.eps"
rep
set terminal qt
