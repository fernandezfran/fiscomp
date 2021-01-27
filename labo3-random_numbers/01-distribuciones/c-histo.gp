#para plotear el histograma de la distribución potencial
set encoding utf8

set xlabel "bin"
set ylabel ""
set grid

array eti[4]
eti[1] = 'ran0'
eti[2] = 'ran2'
eti[3] = 'mzran'
eti[4] = 'M-T'

plot 0.5*exp(-0.5*x) w l, \
  for [i=0:3] 'c-histo.dat' index i u 1:2 w boxes title eti[i+1]

set terminal postscript enhanced color eps
set output "c-histo.eps"
rep
set terminal qt
