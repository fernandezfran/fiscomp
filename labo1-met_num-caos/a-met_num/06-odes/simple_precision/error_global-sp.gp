# grafico en escala log-log del error global en función de 1/h

set xlabel "1/h"
set ylabel "|x(t_f) - x_{exacto}(t_f)|/x_{exacto}(t_f)"
#set title "Problema 6: Error global vs 1/h (log-log)"

set grid
set xrange [1:10**6]
#set yrange [10**(-16):10**(-3)]
set format x "%g"
set format y "10^{%T}"
set log xy

plot 'errores.dat' u (1/($1)):2 smooth unique w lp lt 1 pt 7  lc rgb "red"   t "Euler", \
     ''            u (1/($1)):3 smooth unique w lp lt 1 pt 11 lc rgb "green" t "RK2", \
     ''            u (1/($1)):4 smooth unique w lp lt 1 pt 9  lc rgb "blue"  t "RK4", \
     sqrt(x)*.5*10.**(-7) w l lw 2 dt 2 lc rgb "orange" t "\\~ (1/h)^{1/2}"



set terminal postscript enhanced color eps
set output "error_global-sp.eps"
rep
set terminal qt

