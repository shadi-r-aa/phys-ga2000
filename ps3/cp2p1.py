#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 11:14:52 2023

@author: f003fh5
"""

def f(x):
    return x*(x-1)
def fprime(x):
    return 2*x -1

def deriv(g, x0, delta):
    return 1/delta *(g(x0 + delta) - g(x0))

for i in [10**(-2),10**(-4), 10**(-6), 10**(-8), 10**(-10), 10**(-12), 10**(-14)]:
    print("Derivative is:", deriv(f, 1, i))
    print("Difference is:", abs(fprime(1)-deriv(f, 1, i)))

