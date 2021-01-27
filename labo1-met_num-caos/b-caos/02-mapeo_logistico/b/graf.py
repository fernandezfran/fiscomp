#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

#'tiempo', parte Re e Img de la tranformada de f
t ,retf, imgtf = np.loadtxt('tfx.dat', unpack=True)

t = np.split(t,5)
y = retf * retf + imgtf * imgtf
y = np.split(y,5)
for k in range(0,len(y)):
    y[k] = y[k]/np.max(y[k])


plt.plot(t[4],y[4])
plt.show()
