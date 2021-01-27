#msd vs t
set encoding utf8

set xlabel "t/t_0"
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
#a               = 0.196708         +/- 9.819e-06    (0.004992%)
#b               = 0.0453938        +/- 7.497e-05    (0.1652%)
#
# de donde sale que
# D = 
# en las unidades respectivas de [unidad de long ** 2 / unidad de t]

p 'msdvst.dat' u ($1*9.2424473340241189E-02):2 w p pt 7 t 'simulación', \
  6*x t '\~6t, predicción', \
  0.196708*(x/9.2424473340241189E-02) + 0.0453938 t '\~6D_Lt/D_0, ajuste'

set terminal postscript enhanced color eps 20
set output "difusion.eps"
rep
set terminal qt
