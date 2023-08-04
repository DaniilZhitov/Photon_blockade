import numpy as np

# some default values used as placeholders in simulation functions' arguments.
# many of them (such as OMEGA_DRIVE, U_CONST, LAMBDA_1(2)) should be replaced with values that
# depend on the nature of the simulation
TAU = 2*np.pi
N_CONST = 50
KAPPA_E = TAU*5.0e5
KAPPA_I = TAU*4.0e5
OMEGA_CAVITY = TAU*6e9
U_CONST = TAU*8.0e4
OMEGA_DRIVE = TAU*6e9
LAMBDA_1 = TAU*5e5
LAMBDA_2 = TAU*0