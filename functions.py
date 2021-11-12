import numpy as np
import chaospy as cp

# Get the concentration
def c(x, t, M, n, D, q, Lambda):
    '''Returns the concentration at location x and time t with given porosity (n), dispersion coefficient (D), specific discharge (q) and decay rate (lambda)'''
    c1 = M / np.sqrt(4 * np.pi * D * t)
    c2 = np.exp(-(x - q * t / n) ** 2 / (4 * D * t))
    c3 = np.exp(-Lambda * t)
    return c1 * c2 * c3


# Derivatives for the sensitivity analysis
def dcdn(x, t, M, n, D, q, Lambda):
    return c(x, t, M, n, D, q, Lambda) * q * (- x + q * t / n)/(2 * D * n ** 2)

def dcdD(x, t, M, n, D, q, Lambda):
    return c(x, t, M, n, D, q, Lambda) * (- 1 + (x - q * t / n) ** 2 / (2 * D * t)) / (2  * D)

def dcdq(x, t, M, n, D, q, Lambda):
    return c(x, t, M, n, D, q, Lambda) * (x - q * t / n) / (2 * n * D)

def dcdLambda(x, t, M, n, D, q, Lambda):
    return - t * c(x, t, M, n, D, q, Lambda)  