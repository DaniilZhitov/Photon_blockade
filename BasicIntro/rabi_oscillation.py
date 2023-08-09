import numpy as np
import matplotlib.pyplot as plt
from qutip import *

times = np.linspace(0.0, 20.0, 2000)

psi0 = tensor(fock(2,1))
TAU = 2*np.pi
ratios = np.linspace(0, 5, 25)
n_s = []
for ratio in ratios:
    Gamma = 0.1
    H = 0*sigmaz() + TAU*ratio*Gamma/2 * sigmax()
    result = mesolve(H, psi0, times, [np.sqrt(Gamma)*(sigmax()-1.0j*sigmay()), np.sqrt(.00)*sigmaz()])
    steady_n = expect(result.states[-1], sigmaz())
    n_s.append(steady_n)
plt.ylim(-1,1)
plt.plot(ratios, n_s)
plt.show()
if False:
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