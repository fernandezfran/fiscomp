#script de GNUplot para observar la equilibración de los espines en el 
#    modelos de Ising 2D, en particular, para el inciso a/

set encoding utf8

#set terminal gif animate delay 100
#set output 'MC_eq.gif'

set pm3d map
unset colorb
set size square
set cbrange [-1:1]
set palette defined (-1 'blue', 1 'red')

stats 'configs_eq.dat' u 1 noout
xmax = STATS_max
stats 'configs_eq.dat' u 2 noout
ymax = STATS_max

set xrange [.5 : xmax + .5]
set yrange [.5 : ymax + .5]

steps = STATS_blocks -1

do for [i=0:steps]{
    tt = sprintf("MCS %d", i+1)
    set title tt
    sp 'configs_eq.dat' index i u 1:2:3:3 w image pixels
    pause .2
}
