#trayectoria para distintos factores de crecimiento

set xlabel "t"
set ylabel "x"

set xrange [0:130]
set mxtics 2
set yrange [0:1]
set mytics 4
set key left bottom box opaque
set grid ytics

plot 'xi_0.60-r_1.50.dat' u 1:2 w lp dt 2 pt 5 ps .5 t 'r=1.5',\
     'xi_0.60-r_4.00.dat' u 1:2 w lp dt 2 pt 7 t 'r=4'

set terminal postscript enhanced color eps 26
set termoption enhanced
set output "trayectoria2.eps"
rep
set terminal qt
