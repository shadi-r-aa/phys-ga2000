#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 21:38:34 2023

@author: f003fh5
"""

from random import random
from numpy import arange
from pylab import plot,xlabel,ylabel,show

#Initializing
h=1
tmax = 20000

NBi209 = 0
NPb209 = 0
NTi209 = 0
NBi213 = 10000
pPb209 = 1- 2**(-h/(3.3*60))
pTi209 = 1- 2**(-h/(2.2*60))
pBi213 = 1- 2**(-h/(46*60))




#Decay points and time setting up

tpoints = arange(0.0, tmax, h)
Bi209points = []
Pb209points = []
Ti209points = []
Bi213points = []

#Looping
for t in tpoints:
    Bi209points.append(NBi209)
    Ti209points.append(NTi209)
    Bi213points.append(NBi213)
    Pb209points.append(NPb209)
    
    for i in range(NPb209):
        if random() < pPb209:
            NPb209 -= 1
            NBi209 += 1
    
    for i in range(NTi209):
        if random() < pTi209:
            NTi209 -= 1
            NPb209 += 1
    
    for i in range(NBi213):
        if random() < pBi213:
            NBi213 -= 1
            if random() < 1:
                NPb209 +=1
            else:
                NTi209 +=1

#Plotting

plot(tpoints,Bi209points,label='Bi209')
plot(tpoints,Pb209points,label='Pb209')
plot(tpoints,Ti209points,label='Ti209')
plot(tpoints,Bi213points,label='Bi213')
xlabel('time (seconds)')
ylabel('Number of atoms')
show()


