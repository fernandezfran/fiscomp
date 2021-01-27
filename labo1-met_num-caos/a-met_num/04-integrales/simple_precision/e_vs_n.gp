# grafico en escala log-log el error en función de la cantidad
# de puntos n


set xlabel "N"
set ylabel "{/Symbol e}"
#set title "Problema 4: Error relativo vs N (log-log) - simple precisión"

set grid
set xrange [1:9000]
#set yrange [10**(-16):0.1]
set format x "%g"
set format y "10^{%T}"
set log y
set log x

plot 'errores.dat' u 1:2 w l lt 1 lc rgb "orange" t "Trapecio", \
     '' u 1:3 w l lt 1 lc rgb "green" t "Simpson", \
      sqrt(x)*10.**(-8) w l lw 2 dt 2 lc rgb "blue" t "\\~ N^{1/2}"
#     0.1*x**(-2) w l lw 2 dt 2 lc rgb "orange"  t "\\~ N^{ -2}", \
#     0.01*x**(-4)  w l lw 2 dt 2 lc rgb "green" t "\\~ N^{ -4}"


set terminal postscript enhanced color eps
set output "e_vs_n.eps"
rep
set terminal qt
