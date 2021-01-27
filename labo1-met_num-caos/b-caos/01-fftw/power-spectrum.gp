# grafico en escala de la energía en función del tiempo para los tres métodos

set xlabel "{/Symbol w}" enhanced
set ylabel ""

set grid
set xrange [pi/2.:(1024-1)*2.*pi]  #preguntar?
#set format x "%g"
set format y "%e"
#set key left top               #para mover la leyenda de lugar


plot 'tfx.dat' u (2.*pi/$1):($2**2+$3**2) w l lt 1 lc rgb "red" t "Espectro de potencial"


set terminal postscript enhanced color eps
set termoption enhanced
set output "power-spectrum.eps"
rep
set terminal qt

