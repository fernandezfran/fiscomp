#
set encoding utf8

set xlabel "n pasos"
set ylabel "probabilidad de caer en un cuadrante"
set xrange [-0.1:100.1]
set yrange [.18:.26]
set ytics .18,.02,.26
set mytics 2
set grid
set key r b

p  0.25, \
  "b-ran2-cuadrantes.dat"  u 1:2 w lp pt 11 ps 1.5 dt 2 t "ran2", \
  "b-mzran-cuadrantes.dat" u 1:2 w lp pt 9  ps 1.0 dt 2 t "MZRAN", \
  "b-mt-cuadrantes.dat"    u 1:2 w lp pt 7  ps 0.5 dt 2 t "MT"

set terminal postscript enhanced color eps 18
set output "b-cuadrantes.eps"
rep
set terminal qt
