#gráfico del error en x para t = 180 [s]
set encoding utf8

set xlabel "x [m]" 
set ylabel "|T_{exacta} - T_{numérica}|[{/Symbol \260}C]"
set format y "10^{%T}"
set grid
set log y
set key t c 

plot 'explicito/t1800-explicito-10.dat' u 2:5 w lp dt 2 pt 7 t 'Explícito', \
     'implicito/t1800-implicito-10.dat' u 2:5 w lp dt 2 pt 9 t 'Implícito', \
     'cranck-nicolson/t1800-cn-10.dat'  u 2:5 w lp dt 2 pt 11 t 'Cranck-Nicolson'

set terminal postscript enhanced color eps 26
set output "error_t1800-10.eps"
rep
set terminal qt
