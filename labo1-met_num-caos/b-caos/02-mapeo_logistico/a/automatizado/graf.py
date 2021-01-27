#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def graf(x,y,a,b,title,figname):
    r = ('r = 1.5', 'r = 3.3', 'r = 3.5', 'r = 3.55', 'r = 4')
    plt.title(title)
    plt.xlim(0.0,40.0)
    plt.ylim(-0.1,1.5)
    plt.xlabel('t')
    plt.ylabel('x')
    for i in range(a,b):
        plt.plot(t[i],xt[i],'.--', label=r[i-a])
    plt.legend()
    plt.savefig(figname, dpi=600)
    plt.show()


t, xt = np.loadtxt('trayectorias.dat', unpack=True)

t  = np.split(t, 20)
xt = np.split(xt,20)


graf(t,xt,0,5,'población inicial $x_0$ = 0.06','fig-a-0p06.png')

graf(t,xt,5,10,'población inicial $x_0$ = 0.3','fig-a-0p3.png')

graf(t,xt,10,15,'población inicial $x_0$ = 0.6','fig-a-0p6.png')

graf(t,xt,15,20,'población inicial $x_0$ = 0.9','fig-a-0p9.png')
