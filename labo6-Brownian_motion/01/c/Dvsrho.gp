set encoding utf8

set xlabel '{/Symbol r^{*}}'
set ylabel 'D_L/D_0'
#set xrange [-0.5:12.5]
set grid
set key r t

p 'Dvsrho.dat' u 1:($2/6./9.2424473340241189E-002):($3/6./9.2424473340241189E-002) w errorbar pt 7 t ''

set terminal postscript enhanced color eps 20
set output "Dvsrho.eps"
rep
set terminal qt
