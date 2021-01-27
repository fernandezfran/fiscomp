import numpy as np

#cargo los datos, cada columna en una variable como arreglo
h, e_euler, e_rk2, e_rk4 = np.loadtxt('errores.dat',unpack=True)

#ordeno la primer columna y guardo los indices para poner el resto en este orden
h, idx = np.sort(h), np.argsort(h)

#pongo los errores en el nuevo orden de h
e_euler = e_euler[idx]
e_rk2   = e_rk2[idx]
e_rk4   = e_rk4[idx]

#guardo los datos ordenados en un archivo de datos
a = np.column_stack((h,e_euler,e_rk2, e_rk4))
fmt = '% 1.6e'
np.savetxt('errores-sort.dat', a, fmt=fmt, header='#h, e_euler, e_rk2, e_rk4')
