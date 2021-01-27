#para el calculo del h optimo para la dif. central hay que tener en cuenta que
# \epsilon_a \approx f''(x) * h ** 2 es el error del algoritmo y
# \epsilon_r \approx \frac{\epsilon_maq}{h} es el error de redondeo
# Si sumamos, derivamos e igualamos a cero se obtiene que
# h_opt = (\epsilon_maq/(2. f''(x))) ** (1./3.)
#
# en este caso
#
# f''(1) = e
# \epsilon_maq = ?
#
a = 1.1102230246251565E-016
b = 2. * exp(1)
print "h_opt es igual a: ", (a/b)**(1./3)
