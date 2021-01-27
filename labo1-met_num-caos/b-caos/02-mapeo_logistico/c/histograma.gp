#trayectoria para distintos factores de crecimiento

set xlabel ""
set ylabel "x"

set xrange [0:1]
set xtics 0,.1,1
set mxtics 2 
set grid ytics
set style fill solid 0.25 
set key box opaque at screen 0.9,0.9

plot 'histo.dat' u 1:2 w boxes lc rgb 'green' t 'r=4.0, x_0=0.6'

set terminal postscript enhanced color eps 22
set termoption enhanced
set output "histograma.eps"
rep
set terminal qt
