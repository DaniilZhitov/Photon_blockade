import numpy as np

from CONSTANTS import *
kappa_total = KAPPA_I + KAPPA_E
Lambda_3 = input('Lambda_3/kappa_tot: ')
U = input('U/kappa_i: ')
if Lambda_3 == '':
    Lambda_3 = 0.22*kappa_total
else:
    Lambda_3 = float(Lambda_3)*kappa_total
if U == '':
    U = 0.025*KAPPA_I
else:
    U = float(U)*KAPPA_I
r = 1
alpha = Lambda_3 / (2*U)
Lambda_1 = Lambda_3 * (
                -r + np.abs(Lambda_3) ** 2 / (2 * U ** 2) + 1.0j * kappa_total / (
                   4 * U))
Lambda_2 = -Lambda_3 ** 2 / (4 * U)
Delta = -np.abs(Lambda_3) ** 2 / U

print('Inputs')
print('  L3: '+str(Lambda_3/kappa_total)+'k_tot')
print('  U: '+str(U/KAPPA_I)+'k_i')

print('Step 2 parameters')
print('  alpha: '+'%.2f'%alpha)
print('  L1: '+str(Lambda_1/KAPPA_I)+'k_i')
print('  L2: '+'%.2f'%(Lambda_2/KAPPA_I)+'k_i')
print('  Delta: '+'%.2f'%(Delta/KAPPA_I)+'k_i')

#Limit based on kappa_total * T_drive < 0.02 empirical bound
print('Displacement params')
T_drive = 0.02 / kappa_total
k = alpha / T_drive
Lambda_disp_1_max = 1.0j * k - Delta * alpha
Lambda_disp_2_max = -U*alpha**2
print('  L1_disp_max: '+str(Lambda_disp_1_max/KAPPA_I)+'k_i')
print('  L2_disp_max: '+'%.2f'%(Lambda_disp_2_max/KAPPA_I)+'k_i')


