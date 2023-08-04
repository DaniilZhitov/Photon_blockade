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
        result = BlockadeSimulation.Simulations.blockade_step2(N, U=0.075 * kappa_total, Lambda_3=drive * kappa_total, times=times, dL_1=0.01)
        mean_n = expect(num(N), result.states)
        results.append(result)
    np.savez(
        file='weak_drive_data.npz',
        drives=drives,
        results=results)
def plot_occupation():
    file = np.load('weak_drive_data.npz', allow_pickle=True)
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
    plt.savefig('plot_occupation.png')
def plot_coherence():
    file = np.load('weak_drive_data.npz', allow_pickle=True)
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
    plt.savefig('plot_coherence.png')
#simulate()
#plot_occupation()
plot_coherence()