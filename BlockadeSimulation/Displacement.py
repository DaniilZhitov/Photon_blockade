import time

import numpy as np
from matplotlib import animation
from qutip import *
import matplotlib.pyplot as plt
from CONSTANTS import *

# this file is used for testing different displacement protocols by changing drives, U, Lambda_3
# here we employ a displacement linear in time: alpha = k*t
Lambda_3 = 0.22*KAPPA_I
U = 0.025*KAPPA_I
Delta = -np.abs(Lambda_3) ** 2 / U
N = 120
a = destroy(N)

# displacement rate
k = 1*KAPPA_I
def drive_coeff_1(t, arg):
    return 1.0j*k-Delta*k*t
def drive_coeff_2(t, arg):
    return -U*(k*t)**2


H = [Delta * num(N) + U * a.dag() ** 2 * a ** 2,
     [a.dag(), drive_coeff_1],
     [a, lambda t, arg:drive_coeff_1(t, arg).conjugate()],
     [a.dag()**2, drive_coeff_2],
     [a**2, lambda t, arg:drive_coeff_2(t, arg).conjugate()]]
initial_state = (fock(N, 0)+0.1*fock(N,1)).unit()
times = np.linspace(0, 1/k, 100)
result = mesolve(H, initial_state, times)

X1 = 1/2*expect((a+a.dag()), result.states)
X2 = 1/2*expect(-1j*(a-a.dag()), result.states)
fig, axes = plt.subplots(1, 2, figsize=(12, 3))
x_min, x_max, y_min, y_max = -7, 1, -5, 3

xvec = np.linspace(x_min, x_max, 20)
yvec = np.linspace(y_min, y_max, 20)
flag_plot_wigner = False
if flag_plot_wigner:
    contours = axes[0].contourf(xvec, yvec, wigner(initial_state, xvec, yvec))
    cbar = fig.colorbar(contours)
axes[1].plot(X1, X2)
axes[1].legend()
plt.show()
for ax in axes:
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)


def update(frame):
    print(frame)
    return axes[0].contourf(xvec, yvec, wigner(result.states[frame], xvec, yvec))
flag_animate_wigner = False
if flag_animate_wigner:
    ani = animation.FuncAnimation(fig=fig, func=update, frames=len(result.states), interval=60)
    ani.save(filename='wigner_animation.gif', writer='pillow', fps=15)
#plt.show()