from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from BlockadeSimulation.CONSTANTS import  *
N = 10
U = 0.025*KAPPA_I
a = destroy(N)
x = 0.5*(a + a.dag())
kappa_total = KAPPA_I+KAPPA_E
Lambda_3 = 0.2*kappa_total
Lambda_2 = -Lambda_3 ** 2 / (4 * U)
Delta = np.abs(Lambda_3)**2/(2*U)
H0 = Delta*(num(N)) + U/6 * x ** 3
def drive(t, arg):
    return Lambda_2
initial_state = fock(N, 1)
H = [H0, [a**2+a.dag()**2, drive]]
times = np.linspace(0.0, 1/KAPPA_I, 2000)
result = sesolve(H, initial_state, times)
final_state = result.states[-1]
drive_powers = np.zeros_like(times)
for i in range(len(times)):
    drive_powers[i] = drive(times[i], None)
#plt.plot(times, drive_powers)
plot_fock_distribution(final_state)
plt.show()


