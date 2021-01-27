# grafico del tiempo de ejecución de los dos métodos en función de N

set xlabel "N"
set ylabel "tiempo [s]"

set grid
set key left top


plot 'fftw_vs_dft.dat' u 1:2 w lp pt 7 lt 1 lc rgb "blue"  t "FFTW",\
     ''                u 1:3 w lp pt 7 lt 1 lc rgb "green" t "DFT"


set terminal postscript enhanced color eps
set termoption enhanced
set output "t_vs_N.eps"
rep
set terminal qt

