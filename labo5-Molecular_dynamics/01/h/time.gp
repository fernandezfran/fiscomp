#tiempo vs N
set encoding utf8

set xlabel "N"
set ylabel "t [seg]"
set xrange [100:4100]
set yrange [1:10**3]
set grid
set format x "%g"
set format y "10^{%T}"
set xtics (108,256,500,864,1372,2048, 4000)
set log xy
set key l t

p .125E-3*x**2 t '\~ t^2', \
  '../g/time.dat' u 1:2 w lp pt 7 dt 2 t 'fuerzas O(N^2)', \
  .034*x t '\~ t', \
  'time.dat' u 1:2 w lp pt 5 dt 2 t 'fuerzas O(N)'

set terminal postscript enhanced color eps 20
set output "time.eps"
rep
set terminal qt

