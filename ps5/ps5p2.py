#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 20:46:37 2023

@author: f003fh5
"""

import numpy as np
import pandas as pd
from numpy import exp
import matplotlib.pyplot as plt

def is_float(string):
    """ True if given string is float else False"""
    try:
        return float(string)
    except ValueError:
        return False

data = []
with open('signal.dat', 'r') as f:
    d = f.readlines()
    for i in d:
        k = i.rstrip().split('|')
        for i in k:
            if is_float(i):
                data.append(float(i)) 

data = np.array(data, dtype='float')
time = data[::2]
signal = data[1::2]

plt.scatter(time, signal)
plt.title("signal.dat file")
plt.xlabel("Time (no given units)")
plt.ylabel("Signal (no given units)")
plt.show

#cubic fitting


def cubic_model(x, a, b, c, d):
    return a * x**3 + b * x**2 + c * x + d

A = np.column_stack((time**3, time**2,time, np.ones_like(time)))
U, s, Vt = np.linalg.svd(A, full_matrices=False)
coeff = Vt.T @ np.linalg.inv(np.diag(s)) @ U.T @ signal

a, b, c, d = coeff

signal_cubicfit = cubic_model(time, a, b, c, d)




plt.scatter(time, signal_cubicfit)
plt.title("signal.dat file cubic order fit")
plt.xlabel("Time (no given units)")
plt.ylabel("Signal (no given units)")
plt.show


residuals1= signal - signal_cubicfit
plt.figure(2)
plt.scatter(time, residuals1)
plt.title("signal.dat file cubic order fit residuals")
plt.xlabel("Time (no given units)")
plt.ylabel("Residuals (no given units)")
plt.show

#high polynomial fitting
def polynomial_model(x, *coeffs):
    return np.polyval(coeffs, x)
n=20
A = np.column_stack([time**i for i in range(n, -1, -1)])
U, s, Vt = np.linalg.svd(A, full_matrices=False)
coeff = Vt.T @ np.linalg.inv(np.diag(s)) @ U.T @ signal
signal_20fit = polynomial_model(time, *coeff)

plt.scatter(time, signal_20fit)
plt.title("signal.dat file twentieth order fit")
plt.xlabel("Time (no given units)")
plt.ylabel("Signal (no given units)")
plt.show

residuals2= signal - signal_20fit
plt.figure(2)
plt.scatter(time, residuals2)
plt.title("signal.dat file twentieth order fit residuals")
plt.xlabel("Time (no given units)")
plt.ylabel("Residuals (no given units)")
plt.show


#Harmonic fit

T = (time[-1]-time[0])/2

def fourier_series_with_offset(x,offset, *coeffs):
    n = len(coeffs) // 2
    result = offset
    for i in range(n):
        result += coeffs[i] * np.sin((2 * np.pi / T) * (i+1) * x) + coeffs[i+n] * np.cos((2 * np.pi / T) * (i+1) * x)
    return result

n=40

A = np.column_stack([np.sin((2 * np.pi / T) * (i+1) * time) for i in range(n)] + 
                    [np.cos((2 * np.pi / T) * (i+1) * time) for i in range(n)])

# Perform SVD
U, s, Vt = np.linalg.svd(A, full_matrices=False)

# Solve for the coefficients
coeff = Vt.T @ np.linalg.inv(np.diag(s)) @ U.T @ (signal - np.mean(signal))

# Evaluate the fitted curve
signal_40trigfithalfperiod = fourier_series_with_offset(time, np.mean(signal), *coeff)


plt.figure(3)
plt.scatter(time,signal)
plt.scatter(time, signal_40trigfithalfperiod)
plt.title("signal.dat file 40 fourier terms with period half total elapsed time")
plt.xlabel("Time (no given units)")
plt.ylabel("Signal (no given units)")
plt.show





