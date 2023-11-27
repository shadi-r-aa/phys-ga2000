#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 26 23:23:27 2023

@author: f003fh5
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

L = 1e-8  
x0 = L / 2
sigma = 1e-10
k = 5e10
hbar = 1.0545718e-34

def crank_nicolson(N, T, L, mass, dt, psi0=None):
    dx = L / (N - 1)
    x = np.linspace(0, L, N)
    t = np.linspace(0, T, int(T/dt))

    A = np.empty((3, N-2), dtype=complex)
    B = np.empty((3, N-2), dtype=complex)

    a2 = -0.25j * hbar * dt / (mass * dx**2)
    a1 = 1 + 1j * hbar * dt / (2 * mass * dx**2)
    b2 = 0.25j * hbar * dt / (mass * dx**2)
    b1 = 1 - 1j * hbar * dt / (2 * mass * dx**2)

    A[0, :] = a2 * np.ones(N-2)
    A[1, :] = a1 * np.ones(N-2)
    A[2, :] = a2 * np.ones(N-2)

    B[0, :] = b2 * np.ones(N-2)
    B[1, :] = b1 * np.ones(N-2)
    B[2, :] = b2 * np.ones(N-2)

    if psi0 is None:
        psi0 = np.exp(-(x[1:-1]-x0)**2/2/sigma**2) * np.exp(1j * k * x[1:-1])

    psi = np.zeros((int(T/dt), N), dtype=complex)
    psi[0, 1:N-1] = psi0

    for n in range(int(T/dt) - 1):
        rhs = b1 * psi[n, 1:N-1] + b2 * (psi[n, 2:N] + psi[n, 0:N-2])


        psi[n + 1, 1:N-1] = solve_banded((1, 1), A, rhs)

    return psi, x, t

# Parameters
N = 1000
T = 1e-16
L = 1e-8
mass_electron = 9.10938356e-31
dt = 1e-18

psi_final, x, t = crank_nicolson(N, T, L, mass_electron, dt)

# Plot
for i in range(0, int(T/dt), int(T/dt) // 5):
    psi_normalized = np.real(psi_final[i, 1:N-1]) / np.max(np.abs(psi_final))
    psi_normalized = np.concatenate([[0], psi_normalized, [0]])  # Add zeros at the boundaries
    plt.plot(x, psi_normalized, label=f'Time: {t[i]:.2e}')

plt.title('Normalized Real part of the wavefunction at different times')
plt.xlabel('Position (x)')
plt.ylabel('Normalized Real part of wavefunction')
plt.legend()
plt.show()
