set encoding utf8

set xlabel 'T/T_C'
set ylabel '|M|/N'
set grid
unset key 


p '../b/tmxiecv10.dat' u 1:4:(sqrt(400/9900.)*10.*10.*$6*$1) w errorbar pt 11


set terminal postscript enhanced color eps 18
set output "e.eps"
rep
set terminal qt
