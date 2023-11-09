#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 15:26:50 2023

@author: f003fh5
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
def p(x, beta_0, beta_1):
    return 1/(1+np.exp(-(beta_0+beta_1*x)))
data = pd.read_csv('survey.csv')  
xs = data['age'].to_numpy()
ys = data['recognized_it'].to_numpy()
x_sort = np.argsort(xs)
xs = xs[x_sort]
ys = ys[x_sort]


def log_likelihood(beta, xs, ys):
    beta_0 = beta[0]
    beta_1 = beta[1]
    epsilon = 1e-16
    l_list = [ys[i]*np.log(p(xs[i], beta_0, beta_1)/(1-p(xs[i], beta_0, beta_1)+epsilon)) 
              + np.log(1-p(xs[i], beta_0, beta_1)+epsilon) for i in range(len(xs))]
    ll = np.sum(np.array(l_list), axis = -1)
    return -ll # return log likelihood


from scipy import optimize


xdata = xs
ydata = ys

#Fits
p_start = [0, 0]

# Optimize negative log-likelihood
result = optimize.minimize(log_likelihood, p_start, args=(xs, ys))

# Hessian inverse and residual variance
hess_inv = result.hess_inv
var = result.fun / (len(ys) - len(p_start))

# Covariance matrix of parameters
def Covariance(hess_inv, resVariance):
    return hess_inv * resVariance

# Error of parameters
def error(hess_inv, resVariance):
    covariance = Covariance(hess_inv, resVariance)
    return np.sqrt(np.diag(covariance))

# Calculate errors and covariance matrix
dFit = error(hess_inv, var)
C = Covariance(hess_inv, var)

print('Optimal parameters and error:\n\tp: ', result.x, '\n\tdp: ', dFit)
print('Covariance matrix of optimal parameters:\n\tC: ', C)


b0 = result.x[1]
b1 = result.x[0]

plt.scatter(xs, ys, color='blue', label="Data")

x_values = np.linspace(min(xs), max(xs), 100)

def prob(x):
    return 1 / (1 + np.exp(-x))
 
predicted_probabilities =  prob(b1 +b0*x_values)
plt.plot(x_values, predicted_probabilities, color='red', label='Logistic Regression Model')

plt.xlabel('Age')
plt.ylabel('Probability')
plt.title('Logistic Regression Model')
plt.legend()
plt.show()