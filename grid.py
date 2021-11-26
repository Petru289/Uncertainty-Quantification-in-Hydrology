import numpy as np
import pandas as pd
import chaospy as cp
import matplotlib.pyplot as plt
from functions import c
import parameters as par

N = 10 #number of discretization points on each axis
n = np.linspace(par.nRange[0], par.nRange[1], N)
q = np.linspace(par.qRange[0], par.qRange[1], N)
D = np.linspace(par.DRange[0], par.DRange[1], N)
Lambda = np.linspace(par.LambdaRange[0], par.LambdaRange[1], N)
X, Y, Z, W =  np.meshgrid(n, D, q, Lambda)

#run 'pip install openpyxl' before executing the following line
file = pd.read_excel(r'C:\Users\Gheorghe Pascu\OneDrive - tum.de\WiSe 21-22\Uncertainty Quantification in Hydrology\measurements.xlsx')
obs = file.to_numpy()

t = np.array([5, 10, 15, 20, 25, 30, 35])



vect_c = np.vectorize(c)


def const_on_grid(const, size):
    x = np.full(size, const)
    return np.meshgrid(x, x, x, x)



squares = 0
for i in range(obs.shape[0]):
    squares += (np.full((N, N, N, N), obs[i,0]) - c(5, t[i], par.M, X, Y, Z, W)) ** 2

np.min(squares)
indices = np.where(squares == np.min(squares))

n_opt = n[int(indices[0])]
q_opt = q[int(indices[1])]
D_opt = D[int(indices[2])]
Lambda_opt = Lambda[int(indices[3])]



print(2)






