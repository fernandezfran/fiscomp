set encoding utf8

set xlabel "M"
set ylabel "{/Symbol s}"
#set xrange [-0.5:12.5]
#set yrange [1.0e-7:1.3e-6]
#set ytics 1.0e-7,0.2e-6,1.2e-6 
set grid ytics
set format y "%.1t*10^{%T}"
set key box l t

#set log y
#p '0.1-press' u 1:4:5 w errorbar lc 1 pt 5  t '{/Symbol r} = 0.1', '' u 1:4 w lp noti dt 2 lc 1 pt 5  ,\
#  '0.2-press' u 1:4:5 w errorbar lc 2 pt 7  t '{/Symbol r} = 0.2', '' u 1:4 w lp noti dt 2 lc 2 pt 7  ,\
#  '0.3-press' u 1:4:5 w errorbar lc 3 pt 9  t '{/Symbol r} = 0.3', '' u 1:4 w lp noti dt 2 lc 3 pt 9  ,\
#  '0.4-press' u 1:4:5 w errorbar lc 4 pt 11 t '{/Symbol r} = 0.4', '' u 1:4 w lp noti dt 2 lc 4 pt 11 ,\
#  '0.5-press' u 1:4:5 w errorbar lc 5 pt 13 t '{/Symbol r} = 0.5', '' u 1:4 w lp noti dt 2 lc 5 pt 13 ,\
#  '0.6-press' u 1:4:5 w errorbar lc 6 pt 15 t '{/Symbol r} = 0.6', '' u 1:4 w lp noti dt 2 lc 6 pt 15 ,\
#  '0.7-press' u 1:4:5 w errorbar lc 7 pt 17 t '{/Symbol r} = 0.7', '' u 1:4 w lp noti dt 2 lc 7 pt 17 ,\


p  '0.3-press' u 1:4:5 w errorbar lc 8 pt 19 t '{/Symbol r} = 0.8', '' u 1:4 w lp noti dt 2 lc 8 pt 19 

set terminal postscript enhanced color eps 20
set output "P-ba.eps"
rep
set terminal qt
