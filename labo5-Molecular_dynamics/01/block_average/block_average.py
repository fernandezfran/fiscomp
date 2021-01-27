#!/usr/bin/env python
#
# Script de python para analziar el error en bloques (Appendix D.3 Block Average
#     en el Frenkel) para alguna cantidad A.
#
# ejecutar como:
# $ python3 block_average.py > ba-A.dat
import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt('thermo.dat', usecols=5)
#                                   ^ 1: ekin, 2: epot, 3: etot, 4: temp, 5: pres
print('#numero de bloques, cantidad de datos, media, varianza, error de la varianza')

N_M   = len(A)
media = np.mean(A)
var   = np.var(A)/(N_M-1)
ervar = np.sqrt(2.*var**2/(N_M-1))
print('%2d %6d %e %e %e' % (0, N_M, media, var, ervar))

A_old = A
for k in range(1,13):                                          #loop de bloques
    A_new = np.zeros(int(len(A_old)/2))    #asigno memoria para el nuevo bloque

    for i in range(0,len(A_new)):                     #asigno valores al bloque
        A_new[i] = 0.5*(A_old[2*i - 1] + A_old[2*i])

    N_M   = len(A_new)                     #defino variables que voy a printear
    media = np.mean(A_new)
    var   = np.var(A_new)/(N_M-1)
    ervar = np.sqrt(2.*var**2/(N_M-1))
    print('%2d %6d %e %e %e' % (k, N_M, media, var, ervar))

    A_old = A_new                              #reasigno bloque viejo al actual
