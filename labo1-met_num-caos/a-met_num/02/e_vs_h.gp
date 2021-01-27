# grafico en escala log-log el error en función del h

set encoding utf8

set xlabel "Tamaño de paso = 10^k"
set ylabel "Error"
set title "Error vs h (log-log)"

set grid
set xrange [10**(-15):1]
set yrange [10**(-12):1]
set xtics (1,10**(-5),10**(-7),10**(-10),10**(-15))
unset mytics
set format x "%T"
set format y "10^{%T}"
set log y
set log x

plot 'difnum.dat' u 2:4 w l lt rgb "blue" notitle

set terminal postscript color
set output "e_vs_h.ps"
rep
set terminal qt
