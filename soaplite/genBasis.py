from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object)
from scipy.special import *
from scipy.linalg import *
from scipy.optimize import *
import numpy as np


def myGamma(l, x):
    return gamma(l)*gammaincc(l, x)


def intAllMat(l, a):
    m = np.zeros((a.shape[0], a.shape[0]))
    m[:, :] = a
    m = m + m.transpose()
    return 0.5*gamma(l + 3.0/2.0)*m**(-l-3.0/2.0)


def intAllSqr(l, a):
    return 0.5*gamma(l + 3.0/2.0)*(2*a)**(-l-3.0/2.0)


def intPartSqr(l, a, x):
    return x**(2*l - 1)/2./(2*a)**2*(2*a*x**2)**(0.5 - l)*(gamma(l + 3.0/2.0) - myGamma(l + 3.0/2.0,2*a*x**2))


def minimizeMe(alpha, l, x):
    return np.abs(intPartSqr(l, alpha, x)/intAllSqr(l, alpha) - .99)


def findAlpha(l, a, alphaSpace):
    alphas = np.zeros(a.shape[0])
    for i, j in enumerate(a):
      initG = alphaSpace[np.argmin(minimizeMe(alphaSpace, l, j))]
      alphas[i] = fmin(minimizeMe, x0=initG, args=(l, j), disp=False)
    return alphas


def getOrthNorm(X):
    x = sqrtm(inv(X))
    return x


def getBasisFunc(rcut, n):
    # These are the values for where the different basis functions should decay
    # to: evenly space between 1Å and rcut.
    a = np.linspace(1, rcut, n)
    threshold = 1e-3

    alphas_full = np.zeros((10, n))
    betas_full = np.zeros((10, n, n))

    for l in range(0, 10):
        # The alphas are calculated so that the GTOs will decay to the set
        # threshold value at their respective cutoffs
        alphas = -np.log(threshold/np.power(a, l))/a**2

        # Get the beta factors that orthonormalize the set
        betas = getOrthNorm(intAllMat(l, alphas))

        alphas_full[l, :] = alphas
        betas_full[l, :, :] = betas

    return alphas_full, betas_full


def getBasisFuncSing(rcut, n):
    a = np.linspace(1, rcut, n)
    alphasFull = np.array([])
    betas = np.array([])
    betasFull = np.array([])
    alphaSpace = np.array([])
    val = 0.00001
    for i in range(1100):
        val = val*1.015
        alphaSpace = np.append(alphaSpace, val)
    for l in range(0, 1):
        alphas = findAlpha(l, a, alphaSpace)
        betas = getOrthNorm(intAllMat(l, alphas))
        alphasFull = np.append(alphasFull, alphas)
        betasFull = np.append(betasFull, betas)
    return alphasFull, betasFull


if __name__ == '__main__':
    alphas, betas = getBasisFunc(2.0000, 10)
    np.savetxt('alphasPy.dat', alphas)
    np.savetxt('betasPy.dat',  betas)
