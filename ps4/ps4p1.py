#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 18:25:54 2023

@author: f003fh5
"""

from gaussxw import gaussxwab
import numpy as np
from numpy import linspace, exp
import matplotlib.pyplot as plt

#setting up constants, temperatures, and gaussian quadrature like in the book




def gq(function,N, b):

    xp, wp = gaussxwab(N,a,b)

    s = 0.0

    for k in range(N):
        s += wp[k]*function(xp[k])

    return s

tvals = linspace(5,500,500)
V = 1e-6 
rho = 6.022e28
Td = 428 
kB = 1.38064852e-23 
Nlist= [10,20,30,40,50,60,70]
a = 0
# defining specific heat integrand 

def cv(x):
    return (9*V*rho*kB*(t/Td)**3)*(x**4*exp(x))/(exp(x)-1)**2


lists = list()
for p in range(7):
    lists.append([])
    

for i in range(7):
    for t in tvals:
        b = Td/t
        lists[i].append(gq(cv,Nlist[i],b))

for i in range(7):
    plt.plot(tvals,lists[i], label="N = "+str(Nlist[i]))
plt.title("Heat capacity")
plt.xlabel(r"Temperature (Kelvin)")
plt.ylabel(r"$C_V$ ( Joules per Kelvin)")
plt.xlim(tvals[0],tvals[-1])
plt.legend()
plt.show()



        