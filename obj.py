import numpy as np
import pandas as pd
import parameters as par
from functions import c


t = np.array([5, 10, 15, 20, 25, 30, 35])
file = pd.read_excel('./measurements.xlsx')
obs = file.to_numpy()

def ls(x, y, z, w):
    error = 0
    for i in range(obs.shape[0]):
        error += (obs[i,0] - c(5, t[i], par.M, x, y, z, w)) ** 2 + \
                    (obs[i,1] - c(50, t[i], par.M, x, y, z, w)) ** 2
    return error

def llogs(x, y, z, w):
    error = 0
    for i in range(obs.shape[0]):
        error += (np.log(obs[i,0] + 1) - np.log(c(5, t[i], par.M, x, y, z, w) + 1)) ** 2 + \
                 (np.log(obs[i,1] + 1) - np.log(c(50, t[i], par.M, x, y, z, w) + 1)) ** 2
    return error

def relerr(x, y, z, w):
    error = 0
    for i in range(obs.shape[0]):
        error += ((obs[i,0] - c(5, t[i], par.M, x, y, z, w)) / obs[i,0]) ** 2 + \
                 ((obs[i,1] - c(50, t[i], par.M, x, y, z, w)) / obs[i,1]) ** 2
    return error
    