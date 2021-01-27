set encoding utf8

set ylabel "P(v)"
set xlabel "v"
set grid ytics
set key box
set yrange [0:.35]
set xrange [-5:5]

#p (sqrt(1./2./2./pi))*exp(-x**2/(2.*2.)) lw 1.5 lc 8 t 'Maxwell-Boltzmann',\
#  'histo01p.dat' u 1:2 pt 5 ps 1.2 t '{/Symbol n} = 0.01', \
#  'histo001p.dat' u 1:2 pt 7 ps 0.8 t '{/Symbol n} = 0.001'

p (sqrt(1./2./2./pi))*exp(-x**2/(2.*2.)) lw 1.5 lc 8 t 'Maxwell-Boltzmann',\
  'histo01.dat' u 1:2 pt 5 ps 1.2 t '{/Symbol n} = 0.01', \
  'histo001.dat' u 1:2 pt 7 ps 0.8 t '{/Symbol n} = 0.001'

set term postscript enhanced color eps 22
set output "veldist.eps"
rep
set terminal qt
