from matplotlib import animation
from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import BlockadeSimulation.Simulations
from BlockadeSimulation.CONSTANTS import *
# NOTE: N of Fock states should be at least as high as alpha^2. For our values with alpha ~ 10, this means N > 100
N = 130
kappa_total = KAPPA_I + KAPPA_E
U = 0.025*KAPPA_I
Lambda_3 = 0.15*kappa_total
dL_1 = 0.01
r = 1
alpha = Lambda_3 / (2*U)
# Boundaries of protocol steps
Delta = -np.abs(Lambda_3) ** 2 / U
omega_drive = (OMEGA_CAVITY + 2*U) - Delta
#T_12 = 2/Delta*np.arcsin(np.abs(Delta*alpha/(2*DISP_MAGNITUDE)))
#shifting_one_photon_drive = DISP_MAGNITUDE*np.exp(1.0j*(np.pi/2-Delta*T_12/2))#1.0j*(100*KAPPA_I)#1.0j*alpha / T_12

DISP_MAGNITUDE = 200*KAPPA_I
shifting_one_photon_drive = 1.0j*DISP_MAGNITUDE
T_drive1 = (1.0j*alpha/shifting_one_photon_drive).real
T_drive2 = T_drive1

T_12 = T_drive1
print(kappa_total*T_12)
T_23 = T_12 + 5 / kappa_total
T_f = T_23 + T_drive2
times = np.linspace(0, T_f, 1001)  # the end time is s.t. step 3 is just completed (it starts at T23 and takes T12 to complete)
print('Alpha: '+str(alpha))
frame_fraction = 1001*T_12/(T_f)
print('Frame2: '+str(frame_fraction))

k = alpha / T_drive1
def Lambda_1(t, arg):
    if t < T_12:
        #return shifting_one_photon_drive
        return 1.0j*k-(Delta+2*U)*k*t
    elif t < T_23:
        return Lambda_3 * (
                -r + np.abs(Lambda_3) ** 2 / (2 * U ** 2) + 1.0j * kappa_total / (
                   4 * U)) - dL_1 * Lambda_3
    else:
        #return -shifting_one_photon_drive#*np.exp(-1.0j*Delta*(t-T_23))
        return 1.0j * (-k) - (Delta+2*U) * (k*T_12 - k*(t-T_23))
def Lambda_2(t, arg):
    if t < T_12:
        #return 0
        return -U*(k*t)**2
    elif t < T_23:
        return -Lambda_3 ** 2 / (4 * U)
    else:
        #return 0
        return -U*(k*T_12-k*(t-T_23))**2

def simulate():

    initial_state = fock(N, 0)
    result = BlockadeSimulation.Simulations.resonator_with_time_dependent_drive(
        Lambda_1, Lambda_2, N, KAPPA_E, KAPPA_I, OMEGA_CAVITY, U, omega_drive=omega_drive, initial_state=initial_state, times=times)
    qsave(result, '200_very_advanced_protocol_data')
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
    result = qload('200_very_advanced_protocol_data')
    a = destroy(N)
    X1 = 1 / 2 * expect((a + a.dag()), result.states)
    X2 = -1 / 2 * expect(1j * (a - a.dag()), result.states)
    plt.plot(X1, X2)
    plt.title('Quadrature expectations evolution through the protocol')
    plt.xlabel('$X_1$')
    plt.ylabel('$X_2$')
    plt.savefig('200_very_advanced_plot_trajectory.png')
    #plt.show()
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
    plt.savefig('200_very_advanced_wigner.png')
   # plt.show()
simulate()

result = qload('200_very_advanced_protocol_data')
import math
state12 = result.states[math.ceil(frame_fraction)]
final_state = result.states[len(result.states)-1]
a = destroy(N)
print('g_2 = ' + str(expect(a.dag()**2 * a**2, final_state) / expect(num(N), final_state)**2))


plot_trajectory()
fig, ax = plt.subplots()
plot_fock_distribution(final_state, fig=fig, ax=ax, title='Final state population')
plt.savefig('200_very_advanced_plot_Fock.png')
plot_state_wigner(state12)
#fig, ax = plt.subplots()
#plt.show()
#fig, ax = qutip.visualization.matrix_histogram(final_state.extract_states([0,1,2,3]))
#plt.show()

#plt.xlim(-0.5,10)
#plt.ylim(0,0.6)
#print(expect(num(N), final_state))
#plt.savefig('trajectory1.png')