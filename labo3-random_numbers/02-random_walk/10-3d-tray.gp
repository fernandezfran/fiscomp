# 10 RWs
set encoding utf8

set xlabel "x"
set ylabel "y"
set zlabel "z"
#set xrange [-59:59]
#set mxtics 4
#set yrange [-49:49]
#set ytics -50,20,50
#set mytics 4
set ztics -5,1,5
set grid

splot for [i=0:9] 'random-walk.xyz' index i u 2:3:4 t ''

set terminal postscript enhanced color eps 22
set output "10-3d-tray.eps"
rep
set terminal qt

