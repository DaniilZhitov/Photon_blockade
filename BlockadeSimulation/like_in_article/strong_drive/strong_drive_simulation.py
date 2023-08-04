from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import BlockadeSimulation.Simulations
from BlockadeSimulation.CONSTANTS import *

N = 100
kappa_total = KAPPA_I + KAPPA_E
times = np.linspace(0, 1000 / kappa_total, 10000)
dLs = np.array([0.1, 0.032, 0.01])

def simulate():
    results = []
    for dL in dLs:
        result = BlockadeSimulation.Simulations.blockade_step2(N, U=0.4 * kappa_total, Lambda_3=2 * kappa_total, times=times, dL_1=dL)
        mean_n = expect(num(N), result.states)
       # plt.plot(kappa_total*times, mean_n, label=r'$\tilde{\Lambda}_3=$'+str(drive)+'$\kappa$')
        results.append(result)
    np.savez(
        file='article_strong_drive_data.npz',
        dLs=dLs,
        results=results)
def plot_occupation():
    file = np.load('article_strong_drive_data.npz', allow_pickle=True)
    drives = file['dLs']
    results = file['results']
    for dL, result in zip(drives, results):
        mean_n = expect(num(N), result.states)
        plt.plot(kappa_total * times, mean_n, label=r'$\delta\lambda_1=$' + str(dL))
    plt.xlabel('Time $\kappa t$')
    plt.ylabel(r'$\langle n \rangle$')
    plt.xscale('log')
    plt.xlim(1e-2, 1e3)
    plt.yscale('log')
    plt.ylim(1e-3, 1e2)
    plt.legend()
    plt.show()
    plt.savefig('plot_occupation.png')
def plot_coherence():
    file = np.load('article_strong_drive_data.npz', allow_pickle=True)
    drives = file['dLs']
    results = file['results']
    a = destroy(N)
    for drive, result in zip(drives, results):
        displaced_frame_states = result.states
        g_2 = expect(a.dag()**2 * a**2, displaced_frame_states) / expect(num(N), displaced_frame_states)**2
        plt.plot(kappa_total * times, g_2, label=r'$\tilde{\Lambda}_3=$' + str(drive) + '$\kappa$')
    plt.xlabel('Time $\kappa t$')
    plt.ylabel(r'$g^{(2)}$')
    plt.xscale('linear')
    plt.xlim(0, 2)
    plt.yscale('log')
    plt.ylim(1e-5, 1e2)
    plt.legend()
   # plt.show()
    plt.savefig('plot_coherence.png')
def plot_state_wigner(state):
    fig, ax = plt.subplots()
    x_min, x_max, y_min, y_max = -5, 5, -5, 5
    xvec = np.linspace(x_min, x_max, 30)
    yvec = np.linspace(y_min, y_max, 30)
    quad_mesh = ax.contourf(xvec, yvec, wigner(state, xvec, yvec, method='laguerre'))
    cbar = plt.colorbar(quad_mesh)
    plt.xlabel('Q')
    plt.ylabel('P')
    plt.savefig('WN_state_at_0.8kt.png')
#simulate()
#plot_coherence()

#plot_coherence()
#plot_occupation()

file = np.load('article_strong_drive_data.npz', allow_pickle=True)
drives = file['dLs']
results = file['results']
for dL, result in zip(drives, results):
    state = result.states[8]
    print(expect(num(N),state))
    plot_state_wigner(state)
    break

    #mean_n = expect(num(N), result.states)
    #plt.plot(kappa_total * times, mean_n, label=r'$\delta\lambda_1=$' + str(dL))
