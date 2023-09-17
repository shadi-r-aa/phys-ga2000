#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 20:35:27 2023

@author: f003fh5
"""

import numpy as np

f = np.float32(100.98763)
fint = f.view(np.int32)




fdec = f.item()
diff = 100.98763-np.float32(fdec)




print("The 32-bit float representation :") 
print('{:032b}'.format(fint))

print("The restored value using the mantissa/exponent formula is:")
print(fdec)

print("The difference between the two is:")
print(diff)

