#exponente de Lyapunov y diagrama de orbitas vs r

set xlabel "r"
set ylabel ""
set xrange [3.4:4.0]
set yrange [-0.5:1.2]
set mxtics 4
set mytics 2
set key left top box opaque

p 'x_vs_r-2.dat'      every 3 u 1:2 w p pt 7 ps .025 lc rgb 'gray' t 'diagrama de orbitas', \
  'lyapunov_vs_r.dat' every 2 u 1:2 w lp dt 2 pt 7 ps .5 lc rgb 'red' t '{/Symbol l}'

set terminal postscript enhanced color eps 26
set termoption enhanced
set output "lyap_y_dorb.eps"
rep
set terminal qt
