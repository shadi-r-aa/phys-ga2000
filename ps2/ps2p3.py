#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 21:40:41 2023

@author: f003fh5
"""

#newman 3.7: mandelbrot


import numpy as np
import matplotlib.pyplot as plt

def C(minimum, maximum, finegrain):

    re = np.linspace(minimum, maximum, int((maximum - minimum) * finegrain))

    im = np.linspace(minimum, maximum, int((maximum - minimum) * finegrain))

    return re[np.newaxis, :] + im[:, np.newaxis] * 1j

def conv(c, n_iterate):
    z = 0 
    for i in range(n_iterate):
        z = z**2 + c
    return abs(z) <= 2

def in_set(c, n_iterate):
    element = conv(c, n_iterate)
    return c[element]

c = C(-2,2, 2000)
elements = in_set(c, n_iterate = 100)
plt.xlabel('Re[z]', fontsize=20)
plt.ylabel('Im[z]', fontsize=20)
plt.scatter(elements.real, elements.imag, color = 'black', s=0.01)
