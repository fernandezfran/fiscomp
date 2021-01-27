# grafico en escala log-log del error global en función de 1/h

set xlabel "1/h"
set ylabel "|x(10) - x_{exacto}(10)|/x_{exacto}(10)"
#set title "Problema 6: Error global vs 1/h (log-log)"

set grid
set xrange [1:10**6]
set yrange [10**(-16):10**(2)]
set format x "%g"
set format y "10^{%T}"
set log xy

plot 'errores-sort.dat' u (1/($1)):2  w lp lt 1 pt 7  lc rgb "red"   t "Euler", \
     ''                 u (1/($1)):3  w lp lt 1 pt 11 lc rgb "green" t "RK2", \
     ''                 u (1/($1)):4  w lp lt 1 pt 9  lc rgb "blue"  t "RK4", \
     10*x**(-1) w l lw 2 dt 2 lc rgb "red"   t "\\~ (1/h)^{ -1}", \
     x**(-2) w l lw 2 dt 2 lc rgb "green" t "\\~ (1/h)^{ -2}", \
     0.1*x**(-4) w l lw 2 dt 2 lc rgb "blue"  t "\\~ (1/h)^{ -4}"

set terminal postscript enhanced color eps
set output "error_global.eps"
rep
set terminal qt

