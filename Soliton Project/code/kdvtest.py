#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 12 12:59:44 2023

@author: f003fh5
"""

import unittest
import numpy as np

from module import KdV

class Testkdv(unittest.TestCase):
    def test_solve(self):
        x1=-50
        xN=50
        t1=0
        tN=200
        dx = 0.4
        dt =0.05
        def analytic_solution(x, t):
            return -1/2*np.cosh(1/2*(x-t))**(-2)
        
        def initialcondition_test(x):
            return -1/2*np.cosh(1/2*x)**(-2)

        numerical_result = KdV(dx = dx, dt = dt, eps = -6, mu = 1, x1=x1 ,xN = xN,t1=t1, tN=tN)

        x_values = np.linspace(x1, xN, åå(xN-x1)/dx + 1)
        t_values = np.linspace(t1, tN, (tN-t1)/dx + 1)

        for i, t in enumerate(t_values):
            for j, x in enumerate(x_values):
                self.assertAlmostEqual(numerical_result[i][j], analytic_solution(x, t), places=3)

if __name__ == '__main__':
    unittest.main()