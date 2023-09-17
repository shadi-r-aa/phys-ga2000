#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 17:15:20 2023

@author: f003fh5
"""

import numpy as np
import timeit as timeit

total_sum1 = 0

# Using 3 for loops
L = 100

for i in range(-L,L+1):
    for j in range(-L,L+1):
        for k in range(-L, L+1):
            if i != 0 and j != 0 and k != 0: 
                total_sum1 += (-1)**(i+j+k) *(i**2 + j**2 + k**2)**(-0.5)

%timeit total_sum1
print("Total Sum:", total_sum1)



#Using arrays
def V(i, j, k):
        return (-1)**(i+j+k) *(i**2 + j**2 + k**2)**(-0.5)

indices = np.indices((L,L,L))

mask = np.logical_not(np.all(indices == 0, axis=0))


V_array = V(indices[0], indices[1], indices[2])[mask]

total_sum2 = np.sum(V_array)


%timeit total_sum2
print("Total Sum:", total_sum2) 


#print("Total Sum:", total_sum2)

