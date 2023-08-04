from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import BlockadeSimulation.Simulations
from BlockadeSimulation.CONSTANTS import *

N = 100
kappa_total = KAPPA_I + KAPPA_E
times = np.linspace(0, 100 / kappa_total, 1000)
drives = np.array([0.5, 0.33, 0.22, 0.15, 0.10])

def simulate():
    results = []
    for drive in drives:
        result = BlockadeSimulation.Simulations.blockade_step2(N, U=0.075 * kappa_total, Lambda_3=drive * kappa_total, times=times, dL_1=0.01, r=2)
        mean_n = expect(num(N), result.states)
        results.append(result)
    np.savez(
        file='two_photon_data.npz',
        drives=drives,
        results=results)
def plot_occupation():
    file = np.load('two_photon_data.npz', allow_pickle=True)
    drives = file['drives']
    results = file['results']
    for drive,result in zip(drives, results):
        mean_n = expect(num(N), result.states)
        plt.plot(kappa_total * times, mean_n, label=r'$\tilde{\Lambda}_3=$' + str(drive) + '$\kappa$')
    plt.xlabel('Time $\kappa t$')
    plt.ylabel(r'$\langle n \rangle$')
    plt.xscale('log')
    plt.xlim(1e-1, 1e2)
    plt.yscale('log')
    plt.ylim(1e-4, 1e0)
    plt.legend()
    plt.savefig('two_photon_plotoccupation.png')
def plot_coherence():
    file = np.load('two_photon_data.npz', allow_pickle=True)
    drives = file['drives']
    results = file['results']
    a = destroy(N)
    for drive, result in zip(drives, results):
        displaced_frame_states = result.states
        g_2 = expect(a.dag()**2 * a**2, displaced_frame_states) / expect(num(N), displaced_frame_states)**2
        plt.plot(kappa_total * times, g_2, label=r'$\tilde{\Lambda}_3=$' + str(drive) + '$\kappa$')
    plt.xlabel('Time $\kappa t$')
    plt.ylabel(r'$g^{(2)}$')
    plt.xscale('log')
    plt.xlim(1e-1, 1e2)
    plt.yscale('log')
    plt.ylim(1e-5, 1e1)
    plt.legend()
  #  plt.show()
    plt.savefig('two_photon_plotcoherence.png')
#simulate()
#plot_occupation()
#plot_coherence()

file = np.load('two_photon_data.npz', allow_pickle=True)
drives = file['drives']
results = file['results']
a = destroy(N)
def plot_state_wigner(state):
    fig, ax = plt.subplots()
    x_min, x_max, y_min, y_max = 0, 10, -5, 5
    xvec = np.linspace(x_min, x_max, 30)
    yvec = np.linspace(y_min, y_max, 30)
    quad_mesh = ax.contourf(xvec, yvec, wigner(state, xvec, yvec, method='laguerre'))
    cbar = plt.colorbar(quad_mesh)
    plt.title('Wigner function of the displaced state right before step 2')
    plt.xlabel('Q')
    plt.ylabel('P')
    plt.show()
for drive, result in zip(drives, results):
    #plot_fock_distribution(result.states[len(result.states)-1])
    plot_state_wigner(result.states[len(result.states)-1])
    break
plt.show()
