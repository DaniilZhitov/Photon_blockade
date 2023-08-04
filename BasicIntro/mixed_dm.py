import numpy as np
from qutip import *
import matplotlib.pyplot as plt
# checking Wigner negativity of a mixed 0 and 1 state.
p_1 = 0.45
coh = 0.47
dm = Qobj([[1-p_1, coh], [coh, p_1]])


xvec = np.linspace(-2, 2, 250)
yvec = np.linspace(-2, 2, 250)
mesh = np.meshgrid(xvec, yvec)
wig = wigner(dm, xvec, yvec)

fig, ax = plt.subplots()
quad_mesh = ax.pcolormesh(xvec, yvec, wig)
cbar = plt.colorbar(quad_mesh)
plt.show()
