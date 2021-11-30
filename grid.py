import numpy as np
import pandas as pd
import parameters as par
from functions import c


def grid(N, method): #N = number of discretization points on each axis
    t = np.array([5, 10, 15, 20, 25, 30, 35])

    #run 'pip install openpyxl' before executing the following line
    file = pd.read_excel(r'C:\Users\Gheorghe Pascu\OneDrive - tum.de\WiSe 21-22\Uncertainty Quantification in Hydrology\measurements.xlsx')
    obs = file.to_numpy()

    n = np.linspace(par.nRange[0], par.nRange[1], N)
    D = np.linspace(par.DRange[0], par.DRange[1], N)
    q = np.linspace(par.qRange[0], par.qRange[1], N)
    Lambda = np.linspace(par.LambdaRange[0], par.LambdaRange[1], N)
    X, Y, Z, W =  np.meshgrid(n, D, q, Lambda)
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
    n_opt = n[int(indices[1])]
    D_opt = D[int(indices[0])]
    q_opt = q[int(indices[2])]
    Lambda_opt = Lambda[int(indices[3])]

    return [np.array([n_opt, D_opt, q_opt, Lambda_opt]), min]




if __name__ == "__main__":
    opt_par = grid(30, 'least squares')



