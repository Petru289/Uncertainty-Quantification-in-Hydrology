import numpy as np
import pandas as pd
import chaospy as cp
import parameters as par
from functions import c
import matplotlib.pyplot as plt
import obj


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
        error = obj.ls(X, Y, Z, W)

    if method == 'least log-squares':
        error = obj.llogs(X, Y, Z, W)

    if method == 'relative error':
        error = obj.relerr(X, Y, Z, W)


    min = np.min(error)
    indices = np.where(error == min)
    #VERY CAREFUL HERE, MESHGRID TRICKS YOU. INDICES ARE CORRECT NOW.
    i1, i2, i3, i4 = int(indices[1]), int(indices[0]), int(indices[2]), int(indices[3])
    n_opt = x[i1]
    D_opt = y[i2]
    q_opt = z[i3]
    Lambda_opt = w[i4]

    #Plotting projections of the objective functions

    if method == 'least squares':
        f = obj.ls(X[:, :, i3, i4], Y[:, :, i3, i4], Z[:, :, i3, i4], W[:, :, i3, i4])

    if method == 'least log-squares':
        f = obj.llogs(X[:, :, i3, i4], Y[:, :, i3, i4], Z[:, :, i3, i4], W[:, :, i3, i4])

    if method == 'relative error':
        f = obj.llogs(X[:, :, i3, i4], Y[:, :, i3, i4], Z[:, :, i3, i4], W[:, :, i3, i4])

    levels = [0.0, 0.3, 4.0, 14.0, 38.0, 60.0, 200, 400, 800, 2000, 4000]
    cpl = plt.contour(X[:, :, i3, i4], Y[:, :, i3, i4], f, levels, colors='black', linestyles='dashed', linewidths=1)
    plt.clabel(cpl, inline=1, fontsize=10)
    cpl = plt.contourf(X[:, :, i3, i4], Y[:, :, i3, i4], f, levels, \
                      colors = ['darkgreen', 'green', 'forestgreen', 'seagreen', 'mediumseagreen', \
                                'springgreen', 'aquamarine', 'turquoise', 'paleturquoise', 'lightcyan', 'lightblue'])
    plt.xlabel('n')
    plt.ylabel('D')
    plt.show()

    return [np.array([n_opt, D_opt, q_opt, Lambda_opt]), min]




if __name__ == "__main__":
    opt_par = montecarlo(40, 'relative error')
