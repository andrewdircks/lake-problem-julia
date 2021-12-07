from math import log
import numpy as np
import matplotlib.pyplot as plt


def plot(mu=0.03, sigma=np.sqrt(10**-5), n=10000):
    mc = np.random.lognormal(mean=log(mu**2/np.sqrt(mu**2+sigma**2)), sigma=np.sqrt(log((sigma**2+mu**2)/mu**2)), size=n)
    plt.hist(mc, bins=50)

plot()
plt.yticks([])
plt.xlabel('Natural P Inflows')
plt.show()