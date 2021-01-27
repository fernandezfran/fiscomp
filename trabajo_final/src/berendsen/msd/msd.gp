set encoding utf8

set ylabel "{/Symbol D} r(t)^2"
set xlabel "t"
set key r b box
set log xy

p 'msdvst00.dat' w lp pt 5  ps .75 dt 1 lw 2 t 'NVE',\
  'msdvst0010.dat' w lp pt 7  ps .75 dt 2 lw 2 t '{/Symbol t} = 10',\
  'msdvst0100.dat' w lp pt 9  ps .75 dt 3 lw 2 t '{/Symbol t} = 100',\
  'msdvst1000.dat' w lp pt 11 ps .75 dt 4 lw 2 t '{/Symbol t} = 1000'

set term postscript enhanced color eps 22
set output "msd.eps"
rep
set terminal qt
