#trayectoria para distintos factores de crecimiento

set xlabel "t"
set ylabel "x"

set xrange [0:40]
set xtics 0,10,40
set mxtics 2
set yrange [0:1]
set mytics 4
set key left bottom box opaque
set grid ytics

plot 'xi_0.60-r_1.50.dat' u 1:2 w lp dt 2 pt 5 t 'r=1.5',\
     'xi_0.60-r_3.30.dat' u 1:2 w lp dt 2 pt 7 t 'r=3.3',\
     'xi_0.60-r_3.50.dat' u 1:2 w lp dt 2 pt 9 t 'r=3.5',\
     'xi_0.60-r_3.55.dat' u 1:2 w lp dt 2 pt 11 t 'r=3.55'

set terminal postscript enhanced color eps 26
set termoption enhanced
set output "trayectoria1.eps"
rep
set terminal qt
