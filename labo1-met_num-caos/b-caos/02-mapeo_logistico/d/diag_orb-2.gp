#diagrama de orbitas en el plano x, r

set xlabel "r"
set ylabel "x"
set xrange [3.4:4]
set yrange [0:1]
set mxtics 4
set mytics 2 

p 'x_vs_r-2.dat' every 3 u 1:2 w p pt 7 ps .05 lc rgb 'blue' notitle

set terminal postscript enhanced color eps 26
set termoption enhanced
set output "diag_orb-2.eps"
rep
set terminal qt
