#para plotear el histograma
set encoding utf8

set xlabel "bin"
set ylabel ""
set grid
set yrange [0:.8]
set grid ytics
set style fill transparent solid 0.7

array eti[3]
eti[1] = 'vx'
eti[2] = 'vy'
eti[3] = 'vz'

# T_0 = 1.1
plot for [i=0:2] 'histo.dat' index i u 1:2 w boxes t eti[i+1], \
     sqrt(1./(2.*pi*1.1))*exp(-x**2/(2.*1.1)) lw 1.5 lc 8 t 'distribución de Maxwell'

set terminal pdfcairo color# 22
set termoption enhanced
set output "histo.pdf"
rep
set terminal qt
