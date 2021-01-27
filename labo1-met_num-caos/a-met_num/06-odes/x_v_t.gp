# grafico de la posición y la velocidad en función del tiempo para los distintos
# métodos

set xlabel "t"
set ylabel ""
#set title "Problema 6: Error relativo vs t (log-log)"

set grid
set xrange [0:10]
set yrange [-sqrt(2.) - 0.1 : sqrt(2.) + 0.1]
set mxtics 4
set format x "%g"
set format y "%g"
#set key left top               #para mover la leyenda de lugar


plot 'resultados-euler.dat' u 1:2 w l lt 1 lc rgb "red" t "Euler x(t)", \
     'resultados-rk2.dat'   u 1:2 w l lt 2 lc rgb "green" t "RK2 x(t)", \
     'resultados-rk4.dat'   u 1:2 w l lt 3 lc rgb "blue" t "RK4 x(t)", \
     'resultados-euler.dat' u 1:3 w l dt 1 lc rgb "orange" t "Euler v(t)", \
     'resultados-rk2.dat'   u 1:3 w l dt 2 lc rgb "yellow" t "RK2 v(t)", \
     'resultados-rk4.dat'   u 1:3 w l dt 3 lc rgb "violet" t "RK4 v(t)"

set terminal postscript enhanced color eps
set output "x_v_t.eps"
rep
set terminal qt
