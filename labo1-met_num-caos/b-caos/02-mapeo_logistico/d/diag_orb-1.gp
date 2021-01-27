#diagrama de orbitas en el plano x, r

set xlabel "r"
set ylabel "x"
set xrange [1:4]
set yrange [0:1]
set mxtics 4
set mytics 2 

p 'x_vs_r-1.dat' u 1:2 w p pt 7 ps .1 lc rgb 'blue' notitle

set terminal postscript enhanced color eps
set termoption enhanced
set output "diag_orb-1.eps"
rep
set terminal qt
