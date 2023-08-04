from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from BlockadeSimulation.CONSTANTS import  *
N = 100
g = 40
a = destroy(N)
x = 0.5*(a + a.dag())
omega_c = 100
omega_p = 2*omega_c
signal_amp = 0
drive_amp = 10.0
H0 = omega_c*(num(N)) + g * x ** 4

#PLEASE do not change this file. It took a lot of time to set it up working.

def drive(t, arg):
    raw = np.real(signal_amp*np.exp(1.0j*omega_c*t)+drive_amp*np.exp(1.0j*omega_p*t))
    adiabatic_scale = 0.3*omega_c
    if t < 30/omega_c:
        mid = 15/omega_c
        return raw * np.exp(adiabatic_scale*(t-mid))/(1+np.exp(adiabatic_scale*(t-mid)))
    elif t < 60/omega_c:
        return raw
    else:
        mid = 70/omega_c
        return raw * (1-np.exp(adiabatic_scale*(t-mid))/(1+np.exp(adiabatic_scale*(t-mid))))
initial_state = fock(N, 1)
H = [H0, [a+a.dag(), drive]]
times = np.linspace(0.0, 100.0/omega_c, 2000)
result = sesolve(H, initial_state, times)
final_state = result.states[-1]

drive_powers = np.zeros_like(times)
for i in range(len(times)):
    drive_powers[i] = drive(times[i], None)
#plt.plot(times, drive_powers)
plot_fock_distribution(final_state)
plt.show()


