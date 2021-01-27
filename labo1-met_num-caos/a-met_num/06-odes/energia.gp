# grafico en escala de la energía en función del tiempo para los tres métodos

set xlabel "t"
set ylabel "E"

set grid
set xrange [0:10]
set yrange [0.999998:1.000012]
set mxtics 4
set mytics 2
set format x "%g"
set format y "%t"
set key left top               #para mover la leyenda de lugar


plot 'resultados-euler.dat'every 10000  u 1:4 w l lt 1 lc rgb "red"   t "Euler", \
     'resultados-rk2.dat'  every 1000   u 1:4 w l lt 1 lc rgb "green" t "RK2", \
     'resultados-rk4.dat'  every 10     u 1:4 w l lt 1 lc rgb "blue"  t "RK4"


set terminal postscript enhanced color eps
set output "energia.eps"
rep
set terminal qt

