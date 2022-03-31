import numpy as np
import pandas as pd
import sys
sys.path.insert(0, './Assignment1')
import parameters as par
from functions import c
import obj
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def grid(N, method, plot=False): #N = number of discretization points on each axis
    t = np.array([5, 10, 15, 20, 25, 30, 35])

    n = np.linspace(par.nRange[0], par.nRange[1], N)
    D = np.linspace(par.DRange[0], par.DRange[1], N)
    q = np.linspace(par.qRange[0], par.qRange[1], N)
    Lambda = np.linspace(par.LambdaRange[0], par.LambdaRange[1], N)
    X, Y, Z, W = np.meshgrid(n, D, q, Lambda)

    if method == 'least squares':
        error = obj.ls(X, Y, Z, W)

    if method == 'least log-squares':
        error = obj.llogs(X, Y, Z, W)

    if method == 'relative error':
        error = obj.relerr(X, Y, Z, W)

    minError = np.min(error)
    indices = np.where(error == minError)
    
    #VERY CAREFUL HERE, MESHGRID TRICKS YOU. INDICES ARE CORRECT NOW.
    i1, i2, i3, i4 = int(indices[1]), int(indices[0]), int(indices[2]), int(indices[3])
    n_opt = n[i1]
    D_opt = D[i2]
    q_opt = q[i3]
    Lambda_opt = Lambda[i4]

    #Plotting projections of the objective functions
    if plot:

        # Get X and Y values. If you want to plot different parameters change them here.
        # And don't forget to adjust the labels for the plot as well.
        xValues = X[:, :, i3, i4]
        yValues = Y[:, :, i3, i4]
        zValues = Z[:, :, i3, i4]
        wValues = W[:, :, i3, i4]
        
        # Evalueate the objective function
        f = obj.getError(xValues, yValues, zValues, wValues, method)

        # Define meaningful levels. This highly depends on the objective function.
        if method == 'least squares':
            vmin = 0.02
            levels = np.linspace(vmin, f.max()**(1/6), 10) ** 6
        
        if method == 'least log-squares':
            vmin = 0.02
            levels = np.linspace(vmin, f.max()**(1/6), 10) ** 6
        
        if method == 'relative error':
            vmin = f.min()
            levels = np.linspace(vmin, f.max(), 10)
            levels = np.linspace(vmin, f.max()**(1/10), 10) ** 10

        #levels = [0.02, 0.3, 4.0, 14.0, 38.0, 60.0, 200, 400, 800, 2000, 4000]

        cpl = plt.contour(xValues, yValues, f, levels, colors='black', linestyles='dashed', linewidths=1)
        plt.clabel(cpl, inline=1, fontsize=10)
        cpl = plt.contourf(xValues, yValues, f, levels, cmap='YlGnBu_r', vmin=vmin/2, vmax=f.max(), norm=LogNorm())
        plt.xlabel('Porosity')
        plt.ylabel('Dispersion Coefficient')
        plt.title('Value of the objective function')
        plt.show()

    return [np.array([n_opt, D_opt, q_opt, Lambda_opt]), minError]


if __name__ == "__main__":
    opt_par = grid(30, 'least squares', plot=True)
    print(opt_par)


