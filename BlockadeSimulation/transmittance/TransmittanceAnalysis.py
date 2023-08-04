from qutip import *
import numpy as np
import matplotlib.pyplot as plt
import BlockadeSimulation.Simulations
from BlockadeSimulation.CONSTANTS import *

# transmittance is not really that important

class TransmittanceMeasurement:
    def __init__(self, Delta, Lambda_1, exp_a):
        self.Delta = Delta
        self.Lambda_1 = Lambda_1
        self.exp_a = exp_a

# Physical parameters
kappa_e = KAPPA_E
kappa_i = KAPPA_I
kappa_total = kappa_i + kappa_e
U = 0.2 * kappa_i
omega_cavity = OMEGA_CAVITY

N = 30
times = np.linspace(0, 100/kappa_total, 200)

Deltas = kappa_total*np.linspace(-3, 3, 300)
drive_powers = kappa_total*np.array([0.1, 0.5, 1.0, 1.5, 3.0])
def simulate():
    measurements = []
    for Lambda_1 in drive_powers:
        print('Current drive: '+str(Lambda_1))
        i = 0
        i_total = len(Deltas)
        for Delta in Deltas:
            print('Delta: '+str(i)+'/'+str(i_total))
            i += 1
            # Simulation parameters
            result = BlockadeSimulation.Simulations.resonator(N=N, Lambda_1=Lambda_1, Lambda_2=0,
                                                              initial_state=fock(N, 0), times=times,
                                                              omega_cavity=omega_cavity, omega_drive=omega_cavity-Delta,
                                                              U=U)
            exp_a = expect(destroy(N), result.states[len(result.states)-1])
            measurements.append(TransmittanceMeasurement(Delta, Lambda_1, exp_a))
    measurements = np.array(measurements)
    np.save('transmittance_data.npy', measurements)

def plot_transmittances():
    data = {}
    measurements = np.load('transmittance_data.npy', allow_pickle=True)
    for m in measurements:
        key = m.Lambda_1
        if key in data.keys():
            data[key].append(m)
        else:
            data[key] = [m]
    for drive in data.keys():
        _deltas = []
        _a = []
        data[drive].sort(key=sort_func)

        print(drive / kappa_total)
        for item in data[drive]:
            _deltas.append(item.Delta)
            _a.append(item.exp_a)
        _deltas = np.array(_deltas)
        _a = np.array(_a)
        t = np.abs(1 - 1.0j*(kappa_e/drive) * _a)
        plt.plot(_deltas, t, label=str(drive / kappa_total) + '$\kappa_{tot}$')
    plt.legend()
    plt.xlabel('Detuning '+r'$\Delta/\kappa$')
    plt.ylabel('|t|')
    plt.savefig('plot_transmittance.png')
    #plt.show()
def plot_a_expectation_vs_time():
    Delta = Deltas[150]
    Lambda_1 = drive_powers[0]
    result = BlockadeSimulation.Simulations.resonator(N=N, Lambda_1=Lambda_1, Lambda_2=0,
                                                      initial_state=fock(N, 0), times=times,
                                                      omega_cavity=omega_cavity, omega_drive=omega_cavity - Delta,
                                                      U=U)
    a_exp = expect(destroy(N), result.states)
    plt.plot(times*kappa_total, np.real(a_exp), label='Re')
    plt.plot(times*kappa_total, np.imag(a_exp), label='Im')
    plt.legend()
    plt.show()
def sort_func(m):
    return m.Delta
simulate()
plot_transmittances()





