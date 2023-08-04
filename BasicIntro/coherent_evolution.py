from qutip import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# animation of motion of a coherent state

N = 20
state = coherent(N, 1.3+0.4j)
H = num(N)+0.5
a = destroy(N)
X_1 = 0.5*(a+a.dag())
X_2 = 0.5j * (a.dag()-a)

times = np.linspace(0.0, 10.0, 2000)
result = sesolve(H, state, times, [X_1, X_2])
#plt.plot(result.expect[0], result.expect[1])

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = 4, 3
from matplotlib.animation import FuncAnimation

r = 1 # radius of circle
def quadrature_expectation(frame):
    frame = int(frame)
    return np.array([result.expect[0][frame], result.expect[1][frame]])

# create a figure with an axes
fig, ax = plt.subplots()
# set the axes limits
ax.axis([-1.5, 1.5, -1.5, 1.5])
# set equal aspect such that the circle is not shown as ellipse
ax.set_aspect("equal")
# create a point in the axes
point, = ax.plot(0, 1, marker="o")

# Updating function, to be repeatedly called by the animation
def update(frame):
    # obtain point coordinates
    x,y = quadrature_expectation(frame)
    # set point's coordinates
    point.set_data([x],[y])
    return point,

# create animation with 10ms interval, which is repeated,
# provide the full circle (0,2pi) as parameters
ani = FuncAnimation(fig, update, interval=10, blit=True, repeat=True,
                    frames=np.linspace(0,2000,20000, endpoint=False))

plt.show()