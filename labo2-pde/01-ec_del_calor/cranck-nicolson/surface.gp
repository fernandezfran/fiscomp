#grafico de superficie T(x,t) vs (x,t)

set xlabel "t [s]"
set ylabel "x [m]"
set zlabel "T [{/Symbol \260}C]"
set cblabel "T [{/Symbol \260}C]"

set pm3d
set cbrange [0:100]
set palette rgb 7,5,15 
set colorbox
set hidden3d
set size square

set key opaque
set key at graph 1.5, graph 0.7


splot [0:5500] 'temp_cn.dat' u 1:2:3 w p pt 7 ps .1 t ""


set terminal postscript enhanced color eps 16
set termoption enhanced
set output "surface.eps"
rep
set terminal qt
