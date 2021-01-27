#tiempo vs N
set encoding utf8

set xlabel "N"
set ylabel "t [seg]"
set xrange [100:3000]
set yrange [1:10**3]
set grid
set format x "%g"
set format y "10^{%T}"
set xtics (108,256,500,864,1372,2048)
set log xy
set key l t

p .125E-3*x**2 t '\~ t^2', \
  'time.dat' u 1:2 w lp pt 7 dt 2 t 'mediciones'

set terminal postscript enhanced color eps 20
set output "time.eps"
rep
set terminal qt
