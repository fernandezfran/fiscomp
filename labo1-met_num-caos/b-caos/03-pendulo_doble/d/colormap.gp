#mapa de color para los tiempos de $flip$ 

set xlabel "{/Symbol q_1}"
set ylabel "{/Symbol q_2}"

set pm3d map
set palette defined (0 "green", 10 "#000F00", 10 "#FF0000", \
                     100 "#310000", 100 "purple",1000 "#54025C",1000 "#0000FF", \
                     10000 "#000B70", 10000 "white",10001 "white")
set cbrange [0:10000]
set size ratio -1 #.5
set xrange [-3:3]
set xtics -3,1,3
set yrange [-3:3]
set ytics -3,1,3
set colorbox
unset colorbox

splot '23/t_flip23.dat' every 2 u 1:2:3 noti, '12/t_flip12.dat' every 2 u 1:2:3 noti, \
      '01/t_flip01.dat' every 2 u 1:2:3 noti, \
      '23/t_flip23.dat' every 2 u (-$1):(-$2):3 noti, '12/t_flip12.dat' every 2 u (-$1):(-$2):3 noti, \
      '01/t_flip01.dat' every 2 u (-$1):(-$2):3 noti

set terminal postscript enhanced color eps 20
set termoption enhanced
set output "colormap.gif"
rep
set terminal qt
