# grafico de MSD vs pasos
set encoding utf8

set xlabel "n pasos"
set ylabel "msd"
set xrange [-0.1:1000.1]
set yrange [-50:1150]
#set yrange 
set grid
set key l t

p x w l lw 2 t "N", \
  'a-msd.dat' every 10 u 1:2 w p pt 7 t "msd"

set terminal postscript enhanced color eps 26
set output "a-msd.eps"
rep
set terminal qt
