set encoding utf8

set xlabel 't'
set ylabel 'Fluctuaciones de energía'
set grid
set key l t 

#fluctuaciones totales
p 'thermo.dat' u 1:(abs($2 - 0.420707E+03)) w lp pt 9 ps 1.5 dt 2 t 'cinética', \
  ''           u 1:(abs($3 + 0.139452E+04)) w lp pt 7 ps 0.75 dt 2 t 'potencial', \
  ''           u 1:(abs($4 + 0.973814E+03)) w lp pt 5 dt 2 t 'total'

#fluctuaciones porcentuales
#p 'thermo.dat' u 1:(abs($2 - 0.420707E+03)/0.420707E+03) w lp pt 9 dt 2 t 'cinética', \
#  ''           u 1:(abs($3 + 0.139452E+04)/0.139452E+04) w lp pt 7 dt 2 t 'potencial', \
#  ''           u 1:(abs($4 + 0.973814E+03)/0.973814E+03) w lp pt 5 dt 2 t 'total'

set terminal postscript enhanced color eps 20
set output "energy.eps"
rep
set terminal qt
