#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 25 21:30:14 2023

@author: f003fh5
"""

import numpy as np
from numpy import exp,sqrt
import matplotlib.pyplot as plt
from math import factorial

def H(n,x):
    if n<0:
        return print("Error")
    elif n==0:
        return 1
    elif n==1:
        return 2*x
    else:
        return 2*x*H(n-1,x) - 2*(n-1)*H(n-2,x)

N = 100
xvals= np.linspace(-4,4,100)
nvals = [0,1,2,3]
psivals = np.zeros((4,N))

def psi(n,x):
    return 1/sqrt( 2**n*factorial(n)*np.sqrt(np.pi))*exp(-x**2/2)*H(n,x)

for n in nvals:
    psivals[n,:] = psi(n,xvals)

for n in nvals:
    plt.plot(xvals, psivals[n,:], label="n = "+str(n))
plt.figure(1)
plt.title("Harmonic oscillator wavefunctions")
plt.xlim(-4,4)
plt.xlabel("Position (meters)")
plt.ylabel("Probability amplitude")
plt.legend()
plt.show

plt.figure(2)
Psivals = np.zeros(N)
xvals = np.linspace(-10,10,N)
n = 30
Psivals = psi(n,xvals)
plt.plot(xvals,Psivals,label="n = "+str(n))
plt.xlabel("Position (meters)")
plt.ylabel("Probability amplitude")
plt.xlim(-10,10)
plt.legend()
plt.show()


from gaussxw import gaussxwab
import scipy as sp

N = 100 # points
n = 5
a = -1



def gq(function, b):

    xp, wp = gaussxwab(N,a,b)

    s = 0.0

    for k in range(N):
        s += wp[k]*function(xp[k])

    return s

def ghq(function):

    xhp, whp = sp.special.h_roots(N, mu=False)

    f = 0.0

    for h in range(N):
        f += whp[h]*function(xhp[h])

    return np.float64(f)

#explain transformation 
def ui(z):
    return ((z/(1-z**2))**2)*(1+z**2)/((1-z**2)**2)*abs( psi(n,(z/(1-z**2))))**2

def ui2(y):
    return (y**2)*abs(1/sqrt( 2**n*factorial(n)*np.sqrt(np.pi))*H(n,y))**2


valgq = gq(ui,1)
valghq = ghq(ui2)

print("The difference using GQ is :",np.sqrt(valgq) - 2.3)
print("The difference using GHQ is :",np.sqrt(valghq) - 2.3)

