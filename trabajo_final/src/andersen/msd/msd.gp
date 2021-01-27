set encoding utf8

set ylabel "{/Symbol D} r(t)^2"
set xlabel "t"
set key l t box
set log xy

p 'msdvst00.dat' w lp pt 5  ps .75 dt 1 lw 2 t '{/Symbol n} = 0',\
  'msdvst05.dat' w lp pt 7  ps .75 dt 2 lw 2 t '{/Symbol n} = 5',\
  'msdvst10.dat' w lp pt 9  ps .75 dt 3 lw 2 t '{/Symbol n} = 10',\
  'msdvst20.dat' w lp pt 11 ps .75 dt 4 lw 2 t '{/Symbol n} = 20'

set term postscript enhanced color eps 22
set output "msd.eps"
rep
set terminal qt
