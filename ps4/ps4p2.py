from gaussxw import gaussxwab
import numpy as np
from numpy import linspace, exp
import matplotlib.pyplot as plt

#setting up constants, temperatures, and gaussian quadrature like in the book

N = 20
m = 1
a=0
b=2


xvals = linspace(0.1,2,N)

def gq(function, b):

    xp, wp = gaussxwab(N,a,b)

    s = 0.0

    for k in range(N):
        s += wp[k]*function(xp[k])

    return s



def integrand(x):
    return np.sqrt(8*m)/(np.sqrt(i**4- x**4 ))

L = list()

for i in xvals:
    L.append(gq(integrand,i))
    
plt.title("Anharmonic oscillator")
plt.xlabel("Amplitude (meters)")
plt.ylabel("Period (seconds)")
plt.plot(xvals,L)
plt.show()