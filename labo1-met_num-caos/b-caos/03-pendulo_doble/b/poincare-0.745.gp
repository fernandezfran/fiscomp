#mapa de Poincaré p2 vs theta2

set xlabel "{/Symbol q_2}"
set ylabel "p_2/(m_1L_1^2)"

set xrange [-0.43:0.43]

plot for [i=1:50] 'poincare-0.745.dat' every 2 index i u ($1 - floor(($1 + pi)/2./pi)*2.*pi):2 w p pt 7 ps .2 notitle


set encoding utf8
set terminal postscript enhanced color eps
set termoption enhanced
set output "poincare-0.745.eps"
rep
set terminal qt

