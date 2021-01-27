#Trayectorias de las dos masas
set xlabel "x [m]"
set ylabel "y [m]"

set yrange [-32:]
set xrange [-32:32]

p 'trayectoria-b.dat' every 1000 u (20.*sin($2)):(-20.*cos($2)) w p pt 7 ps .5 lc rgb 'blue' t 'péndulo 1', \
  '' every 100 u (20.*sin($2) + 10.*sin($4)):(-20.*cos($2) - 10.*cos($4)) w p pt 7 ps .25 lc rgb 'red' t 'péndulo 2'


set encoding utf8
set terminal postscript enhanced color eps 20
set termoption enhanced
set output "trayectoria-b.eps"
rep
set terminal qt
