#para plotear el histograma de la distribución potencial
set encoding utf8

set xlabel "bin"
set ylabel ""
set grid

p (1./sqrt(2*pi))*exp(-x**2/2.) w l lw 2 t 'gaussiana', \
  'b-histo.dat' u 1:2 w boxes t 'M-T'

set terminal postscript enhanced color eps
set output "b-histo.eps"
rep
set terminal qt
