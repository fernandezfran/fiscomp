#histograma de la magnetización
#set title "T/T_C = 1.15"
set xlabel "M/N"
set xrange [-1.1:1.1]
set grid ytics
set style fill transparent solid 0.7 
set key r t

plot 'histoM_2.61-128.dat' u 1:2 w boxes lc 4 t 'L = 128', \
     'histoM_2.61-40.dat' u 1:2 w boxes lc 3 t 'L = 40', \
     'histoM_2.61-20.dat' u 1:2 w boxes lc 2 t 'L = 20', \
     'histoM_2.61-10.dat' u 1:2 w boxes lc 1 t 'L = 10'

#set terminal postscript enhanced color eps# 22
set terminal pdfcairo color# 22
set termoption enhanced
set output "histograma-2.61.pdf"
rep
#set terminal qt

