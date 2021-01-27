#histograma de la magnetización
#set title "T/T_C = 0.8"
set xlabel "M/N"
set xrange [-1.1:1.1]
set grid ytics
set style fill transparent solid 0.7 
set key c t

plot 'histoM_2.00-128.dat' u 1:2 w boxes lc 4 t 'L = 128', \
     'histoM_2.00-40.dat' u 1:2 w boxes lc 3 t 'L = 40', \
     'histoM_2.00-20.dat' u 1:2 w boxes lc 2 t 'L = 20', \
     'histoM_2.00-10.dat' u 1:2 w boxes lc 1 t 'L = 10'

#set terminal postscript enhanced color eps# 22
set terminal pdfcairo color# 22
set termoption enhanced
set output "histograma-2.00.pdf"
rep
#set terminal qt

