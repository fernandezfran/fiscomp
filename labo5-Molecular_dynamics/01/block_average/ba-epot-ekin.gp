set encoding utf8

set xlabel 'M'
set ylabel '{/Symbol s}'
set xrange [-0.5:12.5]
set grid
set key l t

p 'ba-ekin.dat' u 1:4:5 w errorbar pt 5 ps 1.1 t 'Energía cinética', \
  'ba-epot.dat' u 1:4:5 w errorbar pt 7 ps 0.9 t 'Energía potencial'

set terminal postscript enhanced color eps 20
set output "ba-epot-ekin.eps"
rep
set terminal qt
