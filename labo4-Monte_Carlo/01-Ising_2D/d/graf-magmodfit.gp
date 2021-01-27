set encoding utf8

set xlabel '|1 - T/T_C|'
set ylabel '|M|/N'
set grid
set key r b 
set xrange [0:0.2]
#set yrange [0.65:0.9]

#f(x) = a*(1 - x)**b
#fit [0:0.1][0.65:0.9] f(x) '../b/tmxiecv128.dat' u (1-$1):4 via a, b


p '../b/tmxiecv128.dat' u (1-$1):(1.21*(abs(1. - $1))**(1./8)) w l t '~ |1-T/T_C|^{1/8}', \
  '' u (1-$1):(abs($4)) w lp pt 5 dt 2 t 'L=128'


set terminal postscript enhanced color eps 18
set output "graf-magmodfit.eps"
rep
set terminal qt
