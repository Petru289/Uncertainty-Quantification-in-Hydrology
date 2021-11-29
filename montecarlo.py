import numpy as np
import pandas as pd
import chaospy as cp
import parameters as par
from functions import c
import matplotlib.pyplot as plt


def montecarlo(N, method):  #N = sample size
 
    t = np.array([5, 10, 15, 20, 25, 30, 35])

    #run 'pip install openpyxl' before executing the following line
    file = pd.read_excel(r'C:\Users\Gheorghe Pascu\OneDrive - tum.de\WiSe 21-22\Uncertainty Quantification in Hydrology\measurements.xlsx')
    obs = file.to_numpy()

    n = cp.Uniform(par.nRange[0], par.nRange[1])
    D = cp.Uniform(par.DRange[0], par.DRange[1])
    q = cp.Uniform(par.qRange[0], par.qRange[1])
    Lambda = cp.Uniform(par.LambdaRange[0], par.LambdaRange[1])
    jointdist = cp.J(n, D, q, Lambda)
    jointsample = jointdist.sample(size = N)
    x = np.sort(jointsample[0,:])
    y = np.sort(jointsample[1,:])
    z = np.sort(jointsample[2,:])
    w = np.sort(jointsample[3,:])
    #This method is maybe inefficient because of sorting the arrays
    X, Y, Z, W = np.meshgrid(x, y, z, w)
    error = 0

    if method == 'least squares':
        for i in range(obs.shape[0]):
            error += (np.full((N, N, N, N), obs[i,0]) - c(5, t[i], par.M, X, Y, Z, W)) ** 2 + \
                     (np.full((N, N, N, N), obs[i,1]) - c(50, t[i], par.M, X, Y, Z, W)) ** 2

    if method == 'least log-squares':
        for i in range(obs.shape[0]):
            error += (np.log(np.full((N, N, N, N), obs[i,0])) - np.log(c(5, t[i], par.M, X, Y, Z, W))) ** 2 + \
                     (np.log(np.full((N, N, N, N), obs[i,1])) - np.log(c(50, t[i], par.M, X, Y, Z, W))) ** 2

    if method == 'relative error':
        for i in range(obs.shape[0]):
            error += ((np.full((N, N, N, N), obs[i,0]) - c(5, t[i], par.M, X, Y, Z, W)) / np.full((N, N, N, N), obs[i,0])) ** 2 + \
                     ((np.full((N, N, N, N), obs[i,1]) - c(50, t[i], par.M, X, Y, Z, W)) / np.full((N, N, N, N), obs[i,1])) ** 2


    min = np.min(error)
    indices = np.where(error == min)
    #VERY CAREFUL HERE, MESHGRID TRICKS YOU. INDICES ARE CORRECT NOW.
    n_opt = x[int(indices[1])]
    D_opt = y[int(indices[0])]
    q_opt = z[int(indices[2])]
    Lambda_opt = w[int(indices[3])]

    # fig, axes = plt.subplots(1, 1)
    # axes.scatter(jointsample[0,:], jointsample[1,:])
    # plt.show()

    print(np.array([n_opt, D_opt, q_opt, Lambda_opt]))




if __name__ == "__main__":
    opt_par = montecarlo(30, 'least squares')
    print(opt_par)
