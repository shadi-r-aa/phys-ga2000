#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 19:11:47 2023

@author: f003fh5
"""
import math
import numpy as np
from scipy import optimize

def brent(f, a, b, tol=1e-9, max_iter=1000):
    # Initialize values
    c = a
    d = 0
    mflag = True
    fm = f(a)
    fa = fm
    fb = f(b)
    if (fa < 0 and fb < 0) or (fa > 0 and fb > 0):
        raise ValueError("No opposite signs at endpoints.")


    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa

    c = a
    fc = fa

    for i in range(max_iter):
        if abs(b-a) < tol:
            return (a + b) / 2

        if fa != fc and fb != fc:
            # Inverse quadratic interpolation
            s = a*fb*fc / ((fa-fb)*(fa-fc)) + b*fa*fc / ((fb-fa)*(fb-fc)) + c*fa*fb / ((fc-fa)*(fc-fb))
        else:
            # Secant
            s = b - fb * (b-a) / (fb-fa)

        if (s < (3*a+b)/4 or s > b) or (mflag and abs(s-b) >= abs(b-c)/2) or (not mflag and abs(s-b) >= abs(c-d)/2) or (mflag and abs(b-c) < tol) or (not mflag and abs(c-d) < tol):
            # Bisection 
            s = (a+b) / 2
            mflag = True
        else:
            mflag = False

        fs = f(s)
        d = c
        c = b

        if (fa*fs < 0):
            b = s
            fb = fs
        else:
            a = s
            fa = fs

        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa

    raise Exception("Doesn't converge within tolerance.")

    
def g(x):
    return np.exp(x)*(x-0.3)*(x+1.7)

myimplementation = brent(g, 0, 1)
scipyimplementation = optimize.brentq(g, 0, 1)
print("Minimum from my implementation:", myimplementation)
print("Minimum from scipy:", scipyimplementation)
    
