set encoding utf8

set xlabel "{/Symbol r}"
set ylabel "P"
set xrange [0.05:0.95]
#set key box l t

p 'Pvsrho.dat' u 1:2:3 w errorbar pt 7 noti,\
  '' u 1:2 w p pt 5 lc 3 ps 2 noti #t 'T = 0.9, {/Symbol n} = 5'

set terminal postscript enhanced color eps 26
set output "Pvsrho.eps"
rep
set terminal qt

