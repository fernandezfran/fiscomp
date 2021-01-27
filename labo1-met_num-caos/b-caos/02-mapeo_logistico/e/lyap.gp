#exponente de Lyapunov vs r

set xlabel "r"
set ylabel "{/Symbol l}"
set mxtics 4
set mytics 2 

p 'lyapunov_vs_r.dat' every 7 u 1:2 w lp dt 2 lc rgb 'red' pt 7 ps .5 notitle

set terminal postscript enhanced color eps 26
set termoption enhanced
set output "lyap.eps"
rep
set terminal qt
