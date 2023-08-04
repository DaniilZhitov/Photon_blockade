from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import BlockadeSimulation.Simulations
from BlockadeSimulation.CONSTANTS import *

N = 20
kappa_total = KAPPA_I + KAPPA_E
times = np.linspace(0, 100 / kappa_total, 10000)
dL = 0.01
Lambda_3 = 1*kappa_total
U = 0.075*U_CONST
r=1

def simulate():
    initial_state = fock(N, 0)
    a = destroy(N)
    H = Lambda_3 * (a.dag()*(num(N)-r) + (num(N)-r)*a)
    dissipators = [np.sqrt(KAPPA_E) * a, np.sqrt(KAPPA_I) * a]

    result = mesolve(H, initial_state, times, dissipators)
    qsave(result, 'kappa=1_strong_result')
def plot_occupation():
    result = qload('kappa=1_strong_result')
    mean_n = expect(num(N), result.states)
    plt.plot(kappa_total * times, mean_n, label=r'$\delta\lambda_1=$' + str(dL))
    plt.xlabel('Time $\kappa t$')
    plt.ylabel(r'$\langle n \rangle$')
    plt.title(r'$\Lambda_3=$'+str(Lambda_3/kappa_total)+r'$\kappa_{total}$')
    #plt.xscale('log')
    #plt.xlim(1e-2, 1e3)
    plt.xlim(0,10)
    #plt.yscale('log')
   # plt.ylim(1e-3, 1e2)
    plt.legend()
    plt.show()
    #plt.savefig('kappa=1_plot_occupation.png')
def plot_coherence():
    file = np.load('kappa=1_strong_drive_data.npz', allow_pickle=True)
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
    plt.savefig('kappa=1_plot_coherence.png')

simulate()
plot_occupation()

#plot_coherence()
