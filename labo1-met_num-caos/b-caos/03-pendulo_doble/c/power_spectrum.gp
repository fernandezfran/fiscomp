#espectro de potencia

set xlabel "{/Symbol w} [Hz]"
set ylabel "espectro de potencia"

set xrange [0:pi]
set yrange [10**(-7):]
#set format x "%.P"
set format y "10^{%T}"
#set key box opaque at screen 0.8,0.9
#set grid
set log y

plot 'power_spectrum-a.dat' u 1:($2**2 + $3**2) w lp pt 7 ps .5 dt 2 lc rgb 'blue' t 'E=-0.6', \
     'power_spectrum-b.dat' u 1:($2**2 + $3**2) w lp pt 7 ps .5 dt 2 lc rgb 'red' t 'E=0'

set terminal postscript enhanced color eps 20
set termoption enhanced
set output "power_spectrum.eps"
rep
set terminal qt
