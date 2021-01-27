# error
set encoding utf8

set xlabel "N"
set ylabel "|I_{numérica} - I_{exacta}|"
set grid
set key l b
set format y "10^{%T}"
set log xy

p .1*x**(-0.5) t 'N^{ -1/2}', \
  'a-error.dat' u 1:4 t 'Monte Carlo', \
  'b-error.dat' u 1:4 t 'Importance Sampling'

set terminal postscript enhanced color eps 20
set output "error-comparacion.eps"
rep
set terminal qt
