#angulos en función del tiempo
set xlabel "t [seg]"
set ylabel "~{/Symbol q_1}{1.1.}, ~{/Symbol q_2}{1.1.}"

#set xrange [0:450]
#set yrange [-1.5:1.5]

p 'trayectoria-a.dat' every 500 u 1:3 w l lt 1 lc rgb 'blue' t 'péndulo 1', \
  '' every 500 u 1:5 w l lt 1 lc rgb 'red' t 'péndulo 2'


set encoding utf8
set terminal postscript enhanced color eps 20
set termoption enhanced
set output "velocidades-a.eps"
rep
set terminal qt
