# grafico de MSD vs pasos
set encoding utf8

set xlabel "n pasos"
set ylabel "msd"
set xrange [-0.1:1000.1]
#set yrange [-50:1150]
#set yrange 
set grid
set key l t

p .01*x w l lw 2 t "0.01 N", \
  'c-msd.dat' every 10 u 1:2 w p pt 7 t "msd"

set terminal postscript enhanced color eps 26
set output "c-msd.eps"
rep
set terminal qt
