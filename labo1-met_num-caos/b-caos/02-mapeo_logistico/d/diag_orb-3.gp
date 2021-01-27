#diagrama de orbitas en el plano x, r

set xlabel "r"
set ylabel "x"
set xtics axis
set xrange [3.847:3.8568]
set xtics 3.847,0.0045,3.8567
set yrange [0.13:0.185]
set mytics 2 

p 'x_vs_r-3.dat' every 2 u 1:2 w p pt 7 ps .1 lc rgb 'blue' notitle

set terminal postscript enhanced color eps 26
set termoption enhanced
set output "diag_orb-3.eps"
rep
set terminal qt
