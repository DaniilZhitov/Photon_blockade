from matplotlib import animation
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import BlockadeSimulation.Simulations
from BlockadeSimulation.CONSTANTS import *
# NOTE: N of Fock states should be at least as high as alpha^2. For our values with alpha ~ 10, this means N > 100
N = 70
kappa_total = KAPPA_I + KAPPA_E
U = 0.2*KAPPA_I
Lambda_3 = 0.33*kappa_total
dL_1 = 0.01
r = 1
alpha = Lambda_3 / (2*U)
# Boundaries of protocol steps
T_12 = 0.001 / kappa_total
T_23 = T_12 + 5 / kappa_total
times = np.linspace(0, T_23+T_12, 10000)  # the end time is s.t. step 3 is just completed (it starts at T23 and takes T12 to complete)

shifting_one_photon_drive = 1.0j*alpha / T_12
def Lambda_1(t, arg):
    if t < T_12:
        return shifting_one_photon_drive
    elif t < T_23:
        return Lambda_3 * (
                -r + np.abs(Lambda_3) ** 2 / (2 * U ** 2) + 1.0j * kappa_total / (
                   4 * U)) - dL_1 * Lambda_3
    else:
        return -shifting_one_photon_drive
def Lambda_2(t, arg):
    if t < T_12:
        return 0
    elif t < T_23:
        return -Lambda_3 ** 2 / (4 * U)
    else:
        return 0

def simulate():
    omega_drive = (OMEGA_CAVITY + 2*U) + np.abs(Lambda_3) ** 2 / U
    initial_state = fock(N, 0)
    result = BlockadeSimulation.Simulations.resonator_with_time_dependent_drive(
        Lambda_1, Lambda_2, N, KAPPA_E, KAPPA_I, OMEGA_CAVITY, U, omega_drive=omega_drive, initial_state=initial_state, times=times)
    qsave(result, 'protocol_data1')
def plot_occupation():
    file = np.load('strong_drive_data.npz', allow_pickle=True)
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
    file = np.load('strong_drive_data.npz', allow_pickle=True)
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
def plot_trajectory():
    result = qload('protocol_data1')
    a = destroy(N)
    X1 = 1 / 2 * expect((a + a.dag()), result.states)
    X2 = -1 / 2 * expect(1j * (a - a.dag()), result.states)
    plt.plot(X1, X2)
    plt.title('Quadrature expectations evolution through the protocol')
    plt.xlabel('$X_1$')
    plt.ylabel('$X_2$')
    #plt.show()
def plot_state_wigner(state):
    fig, ax = plt.subplots()

    x_min, x_max, y_min, y_max = -5, 5, -5, 5
    xvec = np.linspace(x_min, x_max, 30)
    yvec = np.linspace(y_min, y_max, 30)
    quad_mesh = ax.contourf(xvec, yvec, wigner(state, xvec, yvec, method='laguerre'))
    cbar = plt.colorbar(quad_mesh)
    plt.title('Wigner function of the displaced state right before step 2')
    plt.xlabel('Q')
    plt.ylabel('P')
   # plt.show()
#simulate()

result = qload('protocol_data1')

state12 = result.states[1]
final_state = result.states[len(result.states)-1]



plot_trajectory()
fig, ax = plt.subplots()
plot_fock_distribution(final_state, fig=fig, ax=ax, title='Final state population')
fig, ax = qutip.visualization.matrix_histogram_complex(final_state.extract_states([0,1,2,3]))
plot_state_wigner(final_state)
a = destroy(N)
print('g_2 = ' + str(expect(a.dag()**2 * a**2, final_state) / expect(num(N), final_state)**2))
plt.show()
#plt.xlim(-0.5,10)
#plt.ylim(0,0.6)
#print(expect(num(N), final_state))
#plt.savefig('trajectory1.png')