from qutip import *
from BlockadeSimulation.CONSTANTS import *

# Simulating cavity with given parameters and constant drives
def resonator(N=N_CONST, kappa_e=KAPPA_E, kappa_i=KAPPA_I, omega_cavity=OMEGA_CAVITY, U=U_CONST, omega_drive=OMEGA_DRIVE,
                 Lambda_1=LAMBDA_1, Lambda_2=LAMBDA_2, initial_state=None, times=None):
        Delta = (omega_cavity + 2*U) - omega_drive
        a = destroy(N)
        _drive = Lambda_1 * a.dag() + Lambda_2 * a.dag()**2

        H = Delta * num(N) + U * (a.dag() ** 2) * (a ** 2) + _drive + _drive.dag()
        dissipators = [np.sqrt(kappa_e) * a, np.sqrt(kappa_i) * a]
        result = mesolve(H, initial_state, times, dissipators)
        return result

# Simulating step 2 of the blockade with given parameters and constant drives
def blockade_step2(N = N_CONST, kappa_e=KAPPA_E, kappa_i=KAPPA_I, omega_cavity=OMEGA_CAVITY, U=U_CONST,
                   Lambda_3=1, r=1, dL_1 = 0, initial_state=None, times=None):
    kappa_total = kappa_i+kappa_e
    Lambda_1 = Lambda_3 * (
                -r + np.abs(Lambda_3) ** 2 / (2 * U ** 2) + 1.0j * kappa_total / (
                   4 * U)) - dL_1 * Lambda_3
    Lambda_2 = -Lambda_3 ** 2 / (4 * U)
    displace_operator = displace(N, Lambda_3 / (2*U))
    initial_state = displace_operator * fock(N, 0)
    omega_drive = (omega_cavity + 2*U) + np.abs(Lambda_3) ** 2 / U
    result = resonator(N=N, kappa_e=kappa_e, kappa_i=kappa_i, omega_cavity=omega_cavity, U=U, omega_drive=omega_drive, Lambda_1=Lambda_1, Lambda_2=Lambda_2, initial_state=initial_state, times=times)
    displaced_frame_states = []
    for i in range(len(result.states)):
        displaced_frame_states.append(result.states[i].transform(displace_operator.dag()))
    result.states = displaced_frame_states
    return result

# Simulating step 2 of the blockade with given parameters and constant drives
# The difference is that the initial state is not displaced - we start from vacuum
# The steady state should be the same independently of where we start from
def _blockade_step2_start_undisplaced(N = N_CONST, kappa_e=KAPPA_E, kappa_i=KAPPA_I, omega_cavity=OMEGA_CAVITY, U=U_CONST,
                   Lambda_3=1, r=1, dL_1 = 0, initial_state=None, times=None):
    kappa_total = kappa_i+kappa_e
    Lambda_1 = Lambda_3 * (
                -r + np.abs(Lambda_3) ** 2 / (2 * U ** 2) + 1.0j * kappa_total / (
                   4 * U)) - dL_1 * Lambda_3
    Lambda_2 = -Lambda_3 ** 2 / (4 * U)
    displace_operator = displace(N, Lambda_3 / (2*U))
    initial_state = fock(N, 0)
    omega_drive = (omega_cavity + 2*U) + np.abs(Lambda_3) ** 2 / U
    result = resonator(N=N, kappa_e=kappa_e, kappa_i=kappa_i, omega_cavity=omega_cavity, U=U, omega_drive=omega_drive, Lambda_1=Lambda_1, Lambda_2=Lambda_2, initial_state=initial_state, times=times)
    displaced_frame_states = []
    for i in range(len(result.states)):
        displaced_frame_states.append(result.states[i].transform(displace_operator.dag()))
    result.states = displaced_frame_states
    return result

# Self-explanatory. Lambda_1(2) are now functions instead of numbers
def resonator_with_time_dependent_drive(Lambda_1, Lambda_2, N=N_CONST, kappa_e=KAPPA_E, kappa_i=KAPPA_I, omega_cavity=OMEGA_CAVITY, U=U_CONST, omega_drive=OMEGA_DRIVE,
                 initial_state=None, times=None):
    kappa_total = kappa_i + kappa_e
    Delta = (omega_cavity + 2 * U) - omega_drive

    a = destroy(N)
    H0 = Delta * num(N) + U * (a.dag() ** 2) * (a ** 2)

    H = [H0,
         [a.dag(), Lambda_1],
         [a, lambda t, arg : Lambda_1(t, arg).conjugate()],
         [a.dag()**2, Lambda_2],
         [a**2, lambda t, arg: Lambda_2(t, arg).conjugate()]]
    dissipators = [np.sqrt(kappa_e) * a, np.sqrt(kappa_i) * a]

    result = mesolve(H, initial_state, times, dissipators)
    return result
