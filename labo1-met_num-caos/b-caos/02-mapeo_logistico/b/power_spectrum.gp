#espectro de potencia

set xlabel "{/Symbol w} [Hz]"
set ylabel "espectro de potencia"

set xrange [-0.4:3.5]
set yrange [10**(-6):1]
#set format x "%.P"
set format y "10^{%T}"
set key box opaque at screen 0.8,0.9
set grid
set log y

plot 'xi_0.60-r_4.00.dat' u 1:($2**2 + $3**2) w l dt 1 lc rgb 'gray' t 'r=4', \
     'xi_0.60-r_3.55.dat' u 1:($2**2 + $3**2) w lp dt 2 pt 11 ps 1.75 t 'r=3.55', \
     'xi_0.60-r_3.50.dat' u 1:($2**2 + $3**2) w lp dt 2 pt 9 ps 1.75 t 'r=3.5', \
     'xi_0.60-r_3.30.dat' u 1:($2**2 + $3**2) w lp dt 2 pt 7 ps 1.75 t 'r=3.3', \
     'xi_0.60-r_1.50.dat' u 1:($2**2 + $3**2) w lp dt 2 pt 5 ps 1.75 t 'r=1.5'

set terminal postscript enhanced color eps 22
set termoption enhanced
set output "power_spectrum.eps"
rep
set terminal qt
