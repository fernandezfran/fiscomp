#msd vs t
set encoding utf8

set xlabel "t"
set ylabel "MSD"
set grid
set format x "%g"
set format y "10^{%T}"
set log xy
set key l t

# FIT p/coef D
#
#f(x) = a*x + b
#b=1
#a=1
#FIT_LIMIT = 1e-6
#fit [1:10] f(x) 'msdvst-08-nve.dat' via a, b
#
# resultados
#a               = 0.406794         +/- 3.188e-05    (0.007837%)
#b               = -0.00208781      +/- 0.0001938    (9.285%)
#
# de donde sale que
# D = 0.067799 +/- 0.000005
# en las unidades respectivas de [unidad de long ** 2 / unidad de t]

p 'msdvst-08-nve.dat' u 1:2 w lp pt 5 dt 2 t '{/Symbol r} = 0.8', \
  3.*x**2 t '\~ t^2', \
  0.40679*x - 0.0021 t '0.40679 t - 0.0021'

set terminal postscript enhanced color eps 20
set output "difusion.eps"
rep
set terminal qt
