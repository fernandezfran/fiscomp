#mapa de contorno para la ecuación de calor en una barra 1-d
#mostrando las isotermas 

set xlabel "x [m]"
set ylabel "t [s]"
set cblabel "T [{/Symbol \260}C]"

set pm3d map
set cbrange [0:100]
set palette rgb 7,5,15 
set colorbox
set contour surface
set cntrparam levels incr 0,10,100
set size square

set key opaque
set key at graph 1.55, graph 0.7


splot [][0:5500] 'temp_cn.dat' u 2:1:3 t ""


set terminal postscript enhanced color eps 16
set termoption enhanced
set output "colormap.eps"
rep
set terminal qt
