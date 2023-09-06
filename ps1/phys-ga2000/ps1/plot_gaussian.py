# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt


avg = 0;
std = 3;
var = np.square(std);
x= np.arange(-10,10,0.01);
f=np.exp(-np.square(x-avg)/2*var) /np.sqrt(2*np.pi*var);

plt.plot(x,f);
plt.xlabel('x (unit)')
plt.ylabel('y (unit)')
plt.savefig('gaussian.png', )

