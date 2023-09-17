#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 21:55:08 2023

@author: f003fh5
"""

#newman 4.2: quadratic solver

import numpy as np

#part a: defining quadratic solver

def roots1(a,b,c):
    discriminant = b**2 - 4*a*c
    dterm = np.sqrt(abs(discriminant))
    if a ==0 :
        print("Invalid equation")
    if discriminant >0 :
        return((-b + dterm)/(2 * a), (-b - dterm)/(2 * a))
        
    elif discriminant == 0:
        return(-b/(2*a))
    else:
        return((-b + 1j*dterm)/(2 * a),(-b - 1j*dterm)/(2 * a))


    
#part b: changing root expressions
def roots2(a,b,c):
        discriminant = b**2 - 4*a*c
        dterm = np.sqrt(abs(discriminant))
        if a==0:
            print("Invalid equation")
        if discriminant >0 :
            return(2*c/(-b + dterm), 2*c/(-b-dterm))
            
        elif discriminant == 0:
            return(-2*c/(b))
        else:
            return(2*c/(-b + 1j*dterm), 2*c/(-b-1j*dterm))

        


#part c: accurate root finder, add if statement to interpolate between the two depending on relative size of coefficients

def quadratic(a,b,c): 
    x11 = roots1(a,b,c)[0]
    x12 = roots1(a,b,c)[1]
    x21 = roots2(a,b,c)[0]
    x22 = roots2(a,b,c)[1]
    if np.log(abs(x11)) > 0:
        x1 = x11
    else:
        x1 = x21
    if np.log(abs(x12)) > 0:
        x2 = x12
    else:
        x2 = x22
    return (x1,x2)
        
        
   