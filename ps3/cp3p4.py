#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 22:07:15 2023

@author: f003fh5
"""
import numpy as np
from numpy import sort, arange
from pylab import plot,xlabel,ylabel,show
from numpy.random import random

#Definitions and initialization

N = 1000
tau = 3.053*60
mu = np.log(2) / tau

#Generating random numbers
z = random(N)

#Sorting decay times
x = sort(-1/mu * np.log(1-z))
dec = arange(1,N+1)
sur = N - dec

plot(x, sur)
xlabel('time (seconds)')
ylabel('Number of atoms that survived')
show()