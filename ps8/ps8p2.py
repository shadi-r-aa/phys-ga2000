#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 17:44:54 2023

@author: f003fh5
"""
from numpy import array,arange
from pylab import plot,xlabel,ylabel,show


sigma_ = 10
r_=28
b_=8/3

#Defining the equations
def f(r,t):
    x = r[0]
    y = r[1]
    z = r[2]
    fx = sigma_*(y-x)
    fy = r_*x - y - x*z
    fz = x*y - b_*z
    return array([fx,fy,fz],float)

#Time ranges
a = 0.0
b = 100.0
N = 10000
h = (b-a)/N

tpoints = arange(a,b,h)
xpoints = []
ypoints = []
zpoints = []

r = array([0,1,0],float)

for t in tpoints:
    xpoints.append(r[0])
    ypoints.append(r[1])
    zpoints.append(r[2])				
    k1 = h*f(r,t)
    k2 = h*f(r+0.5*k1,t+0.5*h)
    k3 = h*f(r+0.5*k2,t+0.5*h)
    k4 = h*f(r+k3,t+h)
    r += (k1+2*k2+2*k3+k4)/6

plot(tpoints,ypoints)
xlabel("t")
ylabel('y')
show()

plot(zpoints,xpoints)
xlabel('x')
ylabel('z')
show()