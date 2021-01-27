#para plotear fluctuaciones en función de dt
set encoding utf8

set xlabel "r_{cut}"
set ylabel "Fluctuaciones"
set xrange [0.9:5.1]
set yrange [1E-5:1E-1]
set xtics 1,0.5,5
set grid
set format x "%.1t"
set format y "10^{%T}"
set log y
set key r t

plot 'thermo.dat' u 1:5 w lp pt 5 ps 1.5 dt 2 t '{/Symbol D} E_{tot}/{/Symbol D} E_{kin}', \
     ''           u 1:6 w lp pt 7 ps 1 dt 2 t '{/Symbol D} E_{tot}/{/Symbol D} E_{pot}'


set terminal postscript enhanced color eps 20
set output "fluctuaciones.eps"
rep
set terminal qt
