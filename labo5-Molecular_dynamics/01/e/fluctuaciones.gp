#para plotear fluctuaciones en función de dt
set encoding utf8

set xlabel "dt"
set ylabel "{/Symbol D} E_{tot}/{/Symbol D} E_{kin}"
set xrange [:]
set yrange [.25E-5:0.2]
#set xtics (0.1E-04, 0.2E-4, 0.5E-4, 0.1E-03, 0.2E-3, 0.5E-3, 0.1E-2, 0.2E-2, 0.5E-2, 0.1E-1, 0.2E-1)
set grid
set format x "10^{%T}"
set format y "10^{%T}"
set log xy
unset key

plot 'thermo.dat' u 1:4 w lp pt 7 dt 2


set terminal postscript enhanced color eps 20
set output "fluctuaciones.eps"
rep
set terminal qt
