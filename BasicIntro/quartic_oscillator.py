from qutip import *
import numpy as np
import matplotlib.pyplot as plt

N = 200
g = 0.1
a = destroy(N)
x = 0.5*(a + a.dag())
H = (num(N)) + g * x ** 4

initial_state = fock(N, 10)
times = np.linspace(0.0, 10.0, 2000)
result = sesolve(H, initial_state, times, [num(N), x])

plt.xlim(0, 10)
plt.ylim(0, 20)
plt.plot(times, result.expect[0], label='N')
plt.show()


