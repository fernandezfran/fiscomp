# grafico en escala log-log el error en función de la cantidad
# de puntos n

set xlabel "n"
set ylabel "|error|"
set title "Error vs n"

set log y
set log x

plot 'intnum.dat' u 1:2 w l lt rgb "orange" t "Trapecio", '' u 1:3 w l lt rgb "green" t "Simpson"

#ax = GPVAL_DATA_X_MIN + GPVAL_DATA_X_MIN
#bx = GPVAL_DATA_X_MAX + GPVAL_DATA_X_MAX
#ay = GPVAL_DATA_Y_MIN + GPVAL_DATA_Y_MIN
#by = GPVAL_DATA_Y_MAX + GPVAL_DATA_Y_MAX

#set xr [ax:bx]
#set yr [ay:by]



set terminal postscript color
set output "e_vs_n.ps"
rep
set terminal qt
