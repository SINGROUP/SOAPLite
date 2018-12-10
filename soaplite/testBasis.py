from genBasis import getBasisFuncSing
import matplotlib.pyplot as plt
import numpy as np

nMax = 10
alphas, betas = getBasisFuncSing(10.0, nMax)

x = np.linspace(0.01,20,1000)
gss = np.zeros([nMax,len(x)]);



for i in range(nMax):
    plt.plot(x,np.exp(-alphas[i]*x*x))

plt.show()
