# grafico en escala log-log el error local para los distintos métodos
# en función del tiempo


set xlabel "t"
set ylabel "|x(t) - x_{exacto}(t)|"
#set title "Problema 6: Error relativo vs t (log-log)"

set grid
set xrange [0:10]
set yrange [10**(-16):10**(-3)]
set mxtics 4
set format x "%g"
set format y "10^{%T}"
set log y
set key left top               #para mover la leyenda de lugar


plot 'resultados-euler.dat' u 1:5 w l lt 1 lc rgb "red"   t "Euler", \
     'resultados-rk2.dat'   u 1:5 w l lt 1 lc rgb "green" t "RK2", \
     'resultados-rk4.dat'   u 1:5 w l lt 1 lc rgb "blue"  t "RK4"


set terminal postscript enhanced color eps
set output "error_local.eps"
rep
set terminal qt

