import numpy as np
import matplotlib.pyplot as plt
from qutip import *

# a qutip tutorial

times = np.linspace(0.0, 20.0, 2000)

psi0 = tensor(fock(2,1))
TAU = 2*np.pi
H = 0*sigmaz() + np.sqrt(TAU) * 1 * sigmax()

result = mesolve(H, psi0, times, [np.sqrt(0.1)*(sigmax()-1.0j*sigmay()), np.sqrt(.00)*sigmaz()])
state = result.states[61]
if True:
    qutip.visualization.matrix_histogram_complex(state.extract_states([0, 1]))
    plt.show()
    quit()
plt.figure()
plt.plot(times, expect(result.states, sigmaz()))
plt.ylim(-1,1)
plt.xlabel('Time')
plt.ylabel('Expectation values')
plt.legend("cavity photon number")
plt.show()