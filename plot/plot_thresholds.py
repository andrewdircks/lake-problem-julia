import numpy as np
import matplotlib.pyplot as plt


def plot(xs, b, label, q=2.5):
    def simulate(x):
        return (x**q) / (1 + x**q) - b*x
    
    plt.plot(xs, [simulate(x) for x in xs], label=label)


n = 1000
xs = np.linspace(0, 3, num=n)
plt.plot(xs, np.zeros(n))
plot(xs, 0.42, 'Shallow')
plot(xs, 0.62, 'Deep')


plt.xlabel('Lake P Concentration')
plt.ylabel('Net Inward P Flux (Recycling - Removal)')
plt.legend()
plt.show()
