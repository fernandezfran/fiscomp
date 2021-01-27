# 10 RWs
set encoding utf8

set xlabel "x"
set ylabel "y"
set xrange [-59:59]
set mxtics 4
set yrange [-49:49]
set ytics -50,20,50
set mytics 4
set grid

p for [i=0:9] '10-tray.xy' index i using 1:2 w lp ps .1 dt 2 t ''

set terminal postscript enhanced color eps 18
set output "10-tray.eps"
rep
set terminal qt
