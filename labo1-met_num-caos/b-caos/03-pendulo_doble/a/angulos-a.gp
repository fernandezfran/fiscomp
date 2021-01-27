#angulos en función del tiempo
set xlabel "t [seg]"
set ylabel "{/Symbol q_1, q_2}"

set xrange [0:450]
set yrange [-1.5:1.5]

p 'trayectoria-a.dat' every 1000 u 1:2 w l lt 1 lc rgb 'blue' t 'péndulo 1', \
  '' u 1:4 w l lt 1 lc rgb 'red' t 'péndulo 2'


set encoding utf8
set terminal postscript enhanced color eps 20
set termoption enhanced
set output "angulos-a.eps"
rep
set terminal qt
