from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from BlockadeSimulation.CONSTANTS import  *
N = 10
U = 0.025*KAPPA_I
a = destroy(N)
#x = 0.5*(a + a.dag())
omega_c = OMEGA_CAVITY
omega_p = 2*OMEGA_CAVITY
signal_amp = 0
Lambda_3 = 0.2*(KAPPA_I+KAPPA_E)
Lambda_2 = -Lambda_3 ** 2 / (4 * U)
drive_amp = (Lambda_2/U)*omega_c
Delta = np.abs(Lambda_3)**2/(2*U)
def drive(t, arg):
    raw = 000*KAPPA_I+(signal_amp*np.cos(omega_c*t)+drive_amp*np.cos(omega_p*t))
    adiabatic_scale = 30*KAPPA_I
    if t < .4/KAPPA_I:
        mid = .25/KAPPA_I
        return raw * np.exp(adiabatic_scale*(t-mid))/(1+np.exp(adiabatic_scale*(t-mid)))
    elif t < .6/KAPPA_I:
        return raw
    else:
        mid = .75/KAPPA_I
        return raw * (1-np.exp(adiabatic_scale*(t-mid))/(1+np.exp(adiabatic_scale*(t-mid))))
initial_state = fock(N, 1)


a_t = QobjEvo([a, [qeye(N)/omega_c, 'A * cos(w*t)']], args={"A": drive_amp, "w": omega_p})
x_t = a_t + a_t.dag()
x3 = x_t * x_t * x_t
H = (Delta*(a_t.dag()*a_t) + U/6 * x3)
#args = {'A': drive_amp, 'wp': omega_p}
times = np.linspace(0.0, 1/KAPPA_I, 500000)
result = sesolve(H, initial_state, times)
final_state = result.states[-1]

drive_powers = np.zeros_like(times)
for i in range(len(times)):
    drive_powers[i] = drive(times[i], None)
#plt.plot(times, drive_powers)
plot_fock_distribution(final_state)
plt.show()


