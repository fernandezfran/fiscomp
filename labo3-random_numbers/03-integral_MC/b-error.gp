# error
set encoding utf8

set xlabel "N"
set ylabel "|I_{numérica} - I_{exacta}|"
set grid
set key l b
set format y "10^{%T}"
set log xy

p .01*x**(-0.5) t 'N^{ -1/2}', \
  'b-error.dat' u 1:4 t 'error'

set terminal postscript enhanced color eps
set output "b-error.eps"
rep
set terminal qt
