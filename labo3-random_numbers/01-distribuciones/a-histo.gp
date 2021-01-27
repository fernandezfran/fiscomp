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

plot x**(-2) w l, \
  for [i=0:3] 'a-histo.dat' index i u 1:2 w boxes title eti[i+1]

set terminal postscript enhanced color eps
set output "a-histo.eps"
rep
set terminal qt
