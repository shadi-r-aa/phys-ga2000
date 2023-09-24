#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 19 12:42:43 2023

@author: f003fh5
"""
import random as random
import numpy as np
from numpy import zeros
import time
import matplotlib.pyplot as plt




def loopmultiplication(n,A,B):
    C = zeros([n,n], float)
    for i in range(n):
       for j in range(n):
           for k in range(n):
               C[i,j] += A[i,k] * B[k,j]
    
                

def numpydot(n, A,B):
     C = zeros([n,n], float)
     D = np.dot(A,B)


looptimes = np.zeros(0)
for i in [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200]:
    start = time.time()
    loopmultiplication(i,np.random.random((i,i)),np.random.random((i,i)))
    end = time.time()
    looptimes = np.append(looptimes,((end - start) * 1000.0))
    
dottimes = np.zeros(0)
    
for i in [10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200]:
    start = time.time()
    numpydot(i,np.random.random((i,i)),np.random.random((i,i)))
    end = time.time()
    dottimes = np.append(dottimes,((end - start) * 1000.0))




plt.figure()  # Create a new figure
plt.plot([10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200], looptimes, color='blue')
plt.title('Loop method')
plt.xlabel('N')
plt.ylabel('t (microseconds)')

# Create the second plot
plt.figure()  # Create another new figure
plt.plot([10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200], dottimes, color='red')
plt.title('Dot method')
plt.xlabel('N')
plt.ylabel('t (microseconds)')
