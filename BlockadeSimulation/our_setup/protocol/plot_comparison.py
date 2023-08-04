from qutip import *
import matplotlib.pyplot as plt
N = 130


result_advanced = qload('200_very_advanced_protocol_data')
result_simple = qload('simple_protocol_data')
a = destroy(N)
plt.plot(1 / 2 * expect((a + a.dag()), result_simple.states), -1 / 2 * expect(1j * (a - a.dag()), result_simple.states), label='simple')
plt.plot(1 / 2 * expect((a + a.dag()), result_advanced.states), -1 / 2 * expect(1j * (a - a.dag()), result_advanced.states), label='advanced')

plt.title(r'Two displacement protocols')
plt.xlabel('$X_1$')
plt.ylabel('$X_2$')
plt.legend()
plt.savefig('comparison.png')