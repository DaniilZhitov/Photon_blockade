import numpy as np
from qutip import *
import pickle

result = qload('strong_file')
xvec = np.linspace(-3, 3, 20)
yvec = np.linspace(-3, 3, 20)
neg_states = []
for state in result.states:

    wig = wigner(state, xvec, yvec)
    flag = np.any(np.array(wig) < -.0001)
    if flag:
        neg_states.append(state)
np.save('neg_states', neg_states)
print(len(neg_states))