import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, cm
from mpl_toolkits.mplot3d import Axes3D
from qutip import *

# test for animating Wigner function evolution

N = 30
initial_state = coherent(N, 1 / np.sqrt(2) * (1 + 2j))
xvec = np.linspace(-5, 5, 200)
wig = wigner(initial_state, xvec, xvec)

fig, ax = plt.subplots()
my_quad_mesh = ax.pcolormesh(xvec, xvec, wig)

times = np.linspace(0, 10, 100)
a = destroy(N)
H = num(N)
result = mesolve(H, initial_state, times)
def update(frame):
    my_quad_mesh.set_array(wigner(result.states[frame], xvec, xvec))
    return my_quad_mesh


ani = animation.FuncAnimation(fig=fig, func=update, frames=100, interval=30)
plt.show()
quit()




