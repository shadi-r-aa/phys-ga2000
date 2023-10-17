#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  4 21:46:49 2023

@author: f003fh5
"""

import numpy as np
from numpy import exp

import matplotlib.pyplot as plt


def integrand(a,x):
    return x**(a-1) * np.exp(-x)

xvals= np.linspace(0,5,100)
avals = [2,3,4]
integrandvals = np.zeros((5,100))
for n in avals:
    integrandvals[n,:] = integrand(n,xvals)
    
for n in avals:
    plt.plot(xvals, integrandvals[n,:], label="a = "+str(n))
plt.figure(1)
plt.title("Integrand of gamma function")
plt.xlim(0,5)
plt.xlabel("X (no units)")
plt.ylabel("Integrand")
plt.legend()
plt.show

import scipy.integrate as integrate

def gammavalue(a):
    def integranda(z):
        return np.exp( np.log(- (a-1)*z/(z-1))*(a-1) + (a-1)*z/(z-1)) * (a-1)/(z-1)**2
    return integrate.quad(integranda,0,1 )

alist= [ 3/2,3,6,10]

for a in alist:
    print(gammavalue(a))