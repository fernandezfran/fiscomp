#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np

rho_6_01 = np.loadtxt("thermo.6.001000.dat", unpack=True)
rho_6_10 = np.loadtxt("thermo.6.010000.dat", unpack=True)
rho_3_01 = np.loadtxt("thermo.3.001000.dat", unpack=True)
rho_3_10 = np.loadtxt("thermo.3.010000.dat", unpack=True)

# mean = x.sum()/len(x)
# std = sqrt(mean(abs(x - x.mean())**2))
print("rho_0 = 0.6; tau = 1000")
print("valor medio = %f" % (np.mean(rho_6_01[8][int(0.5*len(rho_6_01[8])):])))
print("desviación estándar = %f" % (np.std(rho_6_01[8][int(0.5*len(rho_6_01[8])):])))
print("")

print("rho_0 = 0.6; tau = 10000")
print("valor medio = %f" % (np.mean(rho_6_10[8][int(0.5*len(rho_6_10[8])):])))
print("desviación estándar = %f" % (np.std(rho_6_10[8][int(0.5*len(rho_6_10[8])):])))
print("")

print("rho_0 = 0.3; tau = 1000")
print("valor medio = %f" % (np.mean(rho_3_01[8][int(0.5*len(rho_3_01[8])):])))
print("desviación estándar = %f" % (np.std(rho_3_01[8][int(0.5*len(rho_3_01[8])):])))
print("")

print("rho_0 = 0.3; tau = 10000")
print("valor medio = %f" % (np.mean(rho_3_10[8][int(0.5*len(rho_3_10[8])):])))
print("desviación estándar = %f" % (np.std(rho_3_10[8][int(0.5*len(rho_3_10[8])):])))
print("")
