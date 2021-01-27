set encoding utf8

set xlabel '|T/T_C-1|'
set ylabel 'Suceptibilidad magnética'
set grid
set key r t 
set xrange [0:.5]
set yrange [0:250] 

#f(x) = a*(x - 1)**b
#fit [0:.5][0:250] f(x) '../b/tmxiecv128.dat' u (1-$1):4 via a, b


p '../b/tmxiecv128.dat' u ($1-1):(.18*(abs($1-1)**(-7./4))) w l t '\~|T/T_C-1|^{-7/4}', \
  '' u ($1-1):(128*128*$6) w lp pt 5 dt 2 t 'L=128'


set terminal postscript enhanced color eps 18
set output "graf-chifit.eps"
rep
set terminal qt
