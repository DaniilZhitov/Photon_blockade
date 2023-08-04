import time

import numpy as np
from matplotlib import animation
from qutip import *
import matplotlib.pyplot as plt
from CONSTANTS import *

# This is a slightly different displacement protocol, based on rotation of alpha in the quadrature space
# (compare drive coeffs with Displacement.py)
# it is intended to work in quasi linear regime - for t << 1/Delta motion along the circle's arc in almost like along a tangent

Lambda_3 = 1*KAPPA_I
U = 0.025*KAPPA_I
Delta = -np.abs(Lambda_3) ** 2 / U
alpha = Lambda_3 / (2*U)
N = 1700
a = destroy(N)
T = 1/(5*Delta)
Lambda_1 = 1.0j*alpha / T
print(np.abs(Lambda_1/KAPPA_I))
def drive_coeff_1(t, arg):
    return Lambda_1 + 1.0j*(Delta+2*U)*(Lambda_1*t)
def drive_coeff_2(t, arg):
    alpha_t = (np.exp(-1.0j*(Delta+2*U)*t)-1)*Lambda_1/(Delta+2*U)
    return -U*alpha_t**2

H = [Delta * num(N) + U * a.dag() ** 2 * a ** 2,
     [a.dag(), drive_coeff_1],
     [a, lambda t, arg:drive_coeff_1(t, arg).conjugate()],
     [a.dag()**2, drive_coeff_2],
     [a**2, lambda t, arg:drive_coeff_2(t, arg).conjugate()]]
times = np.linspace(0, T, 1000)
initial_state_0 = fock(N, 0)
initial_state_1 = fock(N, 1)
result_0 = sesolve(H, initial_state_0, times)
result_1 = sesolve(H, initial_state_1, times)

X1_0 = 1/2*expect((a+a.dag()), result_0.states)
X2_0 = 1/2*expect(-1j*(a-a.dag()), result_0.states)
X1_1 = 1/2*expect((a+a.dag()), result_1.states)
X2_1 = 1/2*expect(-1j*(a-a.dag()), result_1.states)

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
x_min, x_max, y_min, y_max = 3, 15, -5, 3
xvec = np.linspace(x_min, x_max, 20)
yvec = np.linspace(y_min, y_max, 20)
axes[1].set_ylim(-.1, .5)
axes[1].plot(X1_0, X2_0, label='0')
axes[1].plot(X1_1, X2_1, label='1')
axes[1].legend()
plt.show()
for ax in axes:
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
